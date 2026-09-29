#!/usr/bin/env python3
"""Resumable stages: plan, compile, simulate, fit, summarize, calibrate (NUTS).

`simulate --block-slice k/n` and `fit --task-slice k/n` process every n-th
block / fit starting at k, so several workers can share one run directory:
each result lives in its own folder and is skipped once written.  Compile the
models first (`compile`) so concurrent workers never race on an executable.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from study import (HERE, ROOT, POPS, SOURCES, MODEL_NAMES, CANDIDATES, load_config, edges, n_haploid,
                   candidates, selections, tasks, msprime_demography, true_ibd, hapibd, snp_summaries,
                   aggregate, stan_data, save_json, event_parameters)
from methods import fitting, ibd, models, parallel, simulate as sim
from methods.blocks import load_blocks
from methods.diagnostics import fit_summary
from methods.fitting import parameter_summary
from methods.utils import file_hashes, open_run as _open_run


def manifest(c):
    files = [HERE / "study.py", HERE / "run_study.py", *sorted((ROOT / "methods").glob("*.py")),
             ROOT / "methods" / "mrca_scan.cpp", *[spec.path for spec in MODEL_NAMES.values()]]
    return {"config": c, "code_sha256": file_hashes(files, ROOT), "candidates": candidates(),
            "selected_blocks": {str(p): {str(cm): selections(c, p, cm) for cm in c["genome_cm"]}
                                for p in range(c["pools"])}}


def pool_folder(out, pool):
    return out / f"pool_{pool:03d}"


def fit_folder(out, pool, cm, rep, source, model, cand, calibration=False):
    return (pool_folder(out, pool) / ("mcmc" if calibration else "fits") / f"cm_{int(cm):04d}"
            / f"rep_{rep:03d}" / source / model / cand)


def simulate(c, a):
    import msprime
    import tskit
    command = a.hapibd_command or ibd.hapibd_command()
    provenance = {"msprime": msprime.__version__, "tskit": tskit.__version__,
                  "numpy": np.__version__, "hapibd_command": command if "hapibd" in a.sources else None}
    if "hapibd" in a.sources:
        provenance["hapibd_version_output"] = ibd.hapibd_version(command)[:5000]
        for item in command:
            if item.endswith(".jar") and Path(item).is_file():
                provenance["jar_sha256"] = hashlib.sha256(Path(item).read_bytes()).hexdigest()
    msdem = msprime_demography(c)
    for pool in range(c["pools"]):
        folder = pool_folder(a.out, pool)
        folder.mkdir(parents=True, exist_ok=True)
        prov_file = folder / ("provenance_" + "_".join(sorted(a.sources)) + ".json")
        if prov_file.exists() and json.loads(prov_file.read_text()) != provenance:
            raise ValueError("Simulator/caller provenance changed; use a new run directory.")
        save_json(prov_file, provenance)
        for block in range(c["blocks"]):
            if a.block_limit is not None and block >= a.block_limit:
                break
            if not parallel.in_slice(block, a.block_slice):
                continue
            cache = folder / f"block_{block:03d}.npz"
            values = dict(np.load(cache)) if cache.exists() else {}
            if all(source + "_count" in values for source in a.sources):
                continue
            seed = sim.block_seed(c["seed"], pool, block)
            tree_file = folder / f"block_{block:03d}.trees"
            t0 = time.monotonic()
            print(f"pool {pool} block {block+1}/{c['blocks']}: ancestry/SNPs", flush=True)
            if tree_file.exists():
                ts = tskit.load(tree_file)
            else:
                ts = sim.simulate_block(msdem, {p: c["diploid_samples"] for p in POPS}, c["block_cm"],
                                        c["recombination_rate"], c["mutation_rate"], seed,
                                        model=c["ancestry_model"])
                ts.dump(tree_file)
            if "snp_numer" not in values:
                values.update(snp_summaries(ts, c))
            for source in a.sources:
                if source + "_count" in values:
                    continue
                print(f"  {source}: extracting segments ({ts.num_trees} marginal trees)", flush=True)
                if source == "true":
                    count, length = true_ibd(ts, c)
                else:
                    count, length = hapibd(ts, c, folder / f"hap_{block:03d}", command)
                values[source + "_count"], values[source + "_length"] = count, length
                values["seed"] = np.asarray(seed)
                tmp = cache.with_suffix(".tmp.npz")
                np.savez_compressed(tmp, **values)
                tmp.replace(cache)
            print(f"  saved {cache.name}, {int(values['snp_sites'].sum())} SNPs; {time.monotonic()-t0:.1f}s",
                  flush=True)


def summaries(sv, weights, obs, candidate, c):
    """Generic fit summary plus recovery of every mapped true parameter."""
    result = fit_summary(sv, weights, obs, n_haploid(c), edges(c))
    result["semantic_parameters"] = {}
    for item in event_parameters(candidate, c):
        stat = parameter_summary(sv[item["variable"]][:, item["index"]], weights)
        stat.update(truth=item["truth"], variable=item["variable"], index=item["index"])
        result["semantic_parameters"][item["name"]] = stat
    if "effective_N" in sv:
        result["semantic_parameters"]["Ne"] = dict(result["effective_N"], truth=c["haploid_ne"],
                                                   variable="effective_N", index=None)
    return result


def fit_stage(c, a, calibration=False):
    by_id = {x["id"]: x for x in candidates()}
    work = tasks(c, a.sources, a.models, a.replicate_limit)
    if calibration:
        work = [t for t in work if t[4] == "T_true" and t[0] in (a.calibration_cm or [c["genome_cm"][-1]])]
    if a.only_cm:
        work = [t for t in work if t[0] in a.only_cm]
    work = [t for i, t in enumerate(work) if parallel.in_slice(i, a.task_slice)]
    compiled = {m: models.normalized_model(MODEL_NAMES[m], a.out / "models") for m in a.models}
    loaded = {}
    for cm, rep, source, model, cand_id in work:
        pool = 0
        folder = fit_folder(a.out, pool, cm, rep, source, model, cand_id, calibration)
        result_path = folder / "result.json"
        if result_path.exists() and not a.retry_failed:
            continue
        if result_path.exists() and json.loads(result_path.read_text()).get("status") == "ok":
            continue
        if pool not in loaded:
            loaded[pool] = load_blocks(pool_folder(a.out, pool), c["blocks"], sorted({"true", *a.sources}))
        selected = selections(c, pool, cm)[rep]
        obs = aggregate(loaded[pool], selected, source, c)
        cand = by_id[cand_id]; spec = MODEL_NAMES[model]
        data = stan_data(cand, obs, c, model)
        folder.mkdir(parents=True, exist_ok=True)
        print(f"cm={cm} rep={rep} {source} {model} {cand_id}", flush=True)
        if calibration:
            result, sv, weights = fitting.nuts(
                compiled[model], data, [fitting.initial(data, spec, s) for s in (1, 7, 13, 19)],
                c["seed"], folder, warmup=a.warmup, samples=a.samples, parallel_chains=a.parallel_chains)
        else:
            result, sv, weights = fitting.pathfinder(
                compiled[model], data, lambda s: fitting.initial(data, spec, s), c["fit_seeds"], folder,
                paths=c["pathfinder_paths"], draws=c["pathfinder_draws"],
                minimum_ess=c["minimum_importance_ess"], maximum_k=c["maximum_pareto_k"],
                keep_csv=c.get("keep_stan_csv", False))
        result = dict(candidate=cand, model=model, **result)
        if result["status"] == "ok":
            result.update(summaries(sv, weights, obs, cand, c))
        result.update(pool=pool, genome_cm=cm, replicate=rep, source=source, selected_blocks=selected)
        save_json(result_path, result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["plan", "compile", "simulate", "fit", "summarize", "calibrate"])
    p.add_argument("--config", type=Path, default=HERE / "config.json")
    p.add_argument("--out", type=Path, default=HERE / "runs" / "default")
    p.add_argument("--sources", nargs="+", choices=SOURCES, default=SOURCES)
    p.add_argument("--models", nargs="+", choices=list(MODEL_NAMES), default=list(MODEL_NAMES))
    p.add_argument("--hapibd-command", type=lambda s: s.split(), default=None,
                   help="argv prefix, e.g. 'java -jar /path/hap-ibd.jar' (default: $HAPIBD_COMMAND or genetics_env)")
    p.add_argument("--block-limit", type=int, help="Debug only: prepare first N blocks")
    p.add_argument("--replicate-limit", type=int, help="Pilot: fit the first N subsamples at every length")
    p.add_argument("--only-cm", type=int, nargs="+", help="Fit only these genome lengths")
    p.add_argument("--task-slice", help="k/n: this worker fits tasks with index % n == k")
    p.add_argument("--block-slice", help="k/n: this worker simulates blocks with index % n == k")
    p.add_argument("--retry-failed", action="store_true")
    p.add_argument("--calibration-cm", type=int, nargs="+", help="calibrate: genome lengths (default: longest)")
    p.add_argument("--warmup", type=int, default=1000)
    p.add_argument("--samples", type=int, default=1000)
    p.add_argument("--parallel-chains", type=int, default=4)
    a = p.parse_args(); a.out = a.out.resolve()
    c = load_config(a.config)
    if a.stage == "summarize":
        # Reading finished results creates no cache, so only the configuration must
        # match the run; later code changes (e.g. to plotting) do not invalidate fits.
        saved = json.loads((a.out / "manifest.json").read_text())
        if saved["config"] != c:
            raise ValueError("Configuration differs from the run being summarised.")
        from visualize import build_summary
        build_summary(c, a.out, a.replicate_limit)
        return
    m = _open_run(a.out, manifest(c))
    if a.stage == "plan":
        n = len(tasks(c, a.sources, a.models))
        print(f"{len(CANDIDATES)} candidates x {len(a.models)} models; SNP-only fitted once per subsample")
        print(f"{c['blocks']} x {c['block_cm']} cM blocks per pool; genome lengths {c['genome_cm']} cM")
        print(f"{c['replicates']} subsamples per length (blocks without replacement); {n} fits "
              f"(each with {len(c['fit_seeds'])} Pathfinder starts)")
        print(f"Manifest: {a.out / 'manifest.json'} ({len(m['candidates'])} candidates)")
    elif a.stage == "compile":
        # Build everything workers share BEFORE they start, so parallel processes
        # never race on the same executable or the MRCA scanner's .so file.
        ibd.mrca_library()
        for model in a.models:
            print(f"{model}: {models.normalized_model(MODEL_NAMES[model], a.out / 'models').exe_file}")
    elif a.stage == "simulate":
        simulate(c, a)
    elif a.stage in ("fit", "calibrate"):
        fit_stage(c, a, calibration=a.stage == "calibrate")


if __name__ == "__main__":
    main()
