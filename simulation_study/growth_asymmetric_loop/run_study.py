#!/usr/bin/env python3
"""Resumable stages: plan, compile, simulate, fit, summarize, calibrate (MCMC).

`simulate --block-slice k/n` and `fit --candidate-slice k/n` process every n-th
block / candidate starting at k, so several workers can share one run
directory: each result lives in its own folder and is skipped once written.
Compile both models first (`compile`) so concurrent workers never race on the
executable.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from study import (HERE, ROOT, POPS, SCENARIOS, SOURCES, VARIANTS, load_config, edges, n_haploid,
                   candidates, selections, msprime_demography, true_ibd, hapibd, snp_summaries,
                   aggregate, stan_data, save_json, event_parameters, branch_map, reference)
from methods import fitting, ibd, models, parallel, simulate as sim
from methods.blocks import load_blocks
from methods.diagnostics import fit_summary
from methods.fitting import parameter_summary
from methods.utils import file_hashes, open_run as _open_run


def manifest(c):
    files = [HERE / "study.py", HERE / "run_study.py", *sorted((ROOT / "methods").glob("*.py")),
             ROOT / "methods" / "mrca_scan.cpp", *[spec.path for spec in VARIANTS.values()]]
    return {"config": c, "code_sha256": file_hashes(files, ROOT),
            "candidate_orders": candidates(),
            "selected_blocks": {str(p): selections(c, p) for p in range(c["pools"])}}


def open_run(c, path):
    return _open_run(path, manifest(c))


def pool_folder(out, scenario, pool):
    return out / scenario / f"pool_{pool:03d}"


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
                import hashlib
                provenance["jar_sha256"] = hashlib.sha256(Path(item).read_bytes()).hexdigest()
    for scenario in a.scenarios:
        msdem = msprime_demography(c, scenario)
        for pool in range(c["pools"]):
            folder = pool_folder(a.out, scenario, pool)
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
                if all(source+"_count" in values for source in a.sources):
                    continue
                seed = sim.block_seed(c["seed"], SCENARIOS.index(scenario), pool, block)
                tree_file = folder / f"block_{block:03d}.trees"
                t0 = time.monotonic()
                print(f"{scenario} pool {pool} block {block+1}/{c['blocks']}: ancestry/SNPs", flush=True)
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
                    if source+"_count" in values:
                        continue
                    print(f"  {source}: extracting segments ({ts.num_trees} marginal trees)", flush=True)
                    if source == "true":
                        count, length = true_ibd(ts, c)
                    else:
                        count, length = hapibd(ts, c, folder / f"hap_{block:03d}", command)
                    values[source+"_count"], values[source+"_length"] = count, length
                    values["seed"] = np.asarray(seed)
                    tmp = cache.with_suffix(".tmp.npz")
                    np.savez_compressed(tmp, **values)
                    tmp.replace(cache)
                print(f"  saved {cache.name}, {int(values['snp_sites'].sum())} SNPs; {time.monotonic()-t0:.1f}s", flush=True)


def truth_for(candidate, c):
    mapping = event_parameters(candidate, c)
    if not mapping:
        return None
    times = [x["truth"] for x in mapping if x["variable"] == "cumulative_times"]
    return dict(times=times, correct_order=bool(np.all(np.diff(times) > 0)))


def summaries(sv, weights, obs, candidate, c, scenario):
    """Generic fit summary plus this study's semantic parameters and Ne references."""
    result = fit_summary(sv, weights, obs, n_haploid(c), edges(c))
    result["semantic_parameters"] = {}
    for item in event_parameters(candidate, c):
        values = sv[item["variable"]][:, item["index"]]
        if item["fold"]:
            values = np.maximum(values, 1-values)
        stat = parameter_summary(values, weights)
        stat.update(truth=item["truth"], variable=item["variable"], index=item["index"], folded=item["fold"])
        result["semantic_parameters"][item["name"]] = stat
    # True-size references for every node of a correct candidate: the interval the
    # branch spans, its sizes at both ends, and the harmonic/arithmetic means of
    # the true effective trajectory.  Loop-omitting candidates get the COLLAPSED
    # composite for b; the explicit-loop graph gets the generating branches.
    bmap = branch_map(candidate)
    result["ne_reference"] = [dict(node=node, **reference(c, scenario, bmap[node],
                                                          collapsed=not candidate["explicit_loop"]))
                              if node in bmap else dict(node=node, branch=None)
                              for node in candidate["nodes"]]
    return result


def fit_stage(c, a, calibration=False):
    chosen = candidates()
    if calibration:
        chosen = [x for x in chosen if x["correct"] and truth_for(x, c)["correct_order"]]
    if a.only_graphs:
        chosen = [x for x in chosen if x["graph"] in a.only_graphs]
    chosen = [x for i, x in enumerate(chosen) if parallel.in_slice(i, a.candidate_slice)]
    if not chosen:
        raise ValueError("No requested candidate graphs")
    for variant in a.variants:
        spec = VARIANTS[variant]
        model = models.normalized_model(spec, a.out / "models")
        for scenario in a.scenarios:
            for pool in range(c["pools"]):
                pf = pool_folder(a.out, scenario, pool)
                block_list = load_blocks(pf, c["blocks"], a.sources)
                for rep, selected in enumerate(selections(c, pool)):
                    if a.replicate_limit is not None and rep >= a.replicate_limit:
                        break
                    for source in a.sources:
                        obs = aggregate(block_list, selected, source, c)
                        for cand in chosen:
                            folder = pf / ("mcmc" if calibration else "fits") / f"rep_{rep:03d}" / source / variant / cand["id"]
                            result_path = folder / "result.json"
                            if result_path.exists() and not a.retry_failed:
                                continue
                            if result_path.exists() and json.loads(result_path.read_text()).get("reliable"):
                                continue
                            folder.mkdir(parents=True, exist_ok=True)
                            data = stan_data(cand, obs, c, variant)
                            print(f"{scenario} pool={pool} rep={rep} {source} {variant} {cand['id']}", flush=True)
                            if calibration:
                                result, sv, weights = fitting.nuts(
                                    model, data, [fitting.initial(data, spec, s) for s in (1, 7, 13, 19)],
                                    c["seed"], folder, warmup=a.warmup, samples=a.samples,
                                    parallel_chains=a.parallel_chains)
                            else:
                                result, sv, weights = fitting.pathfinder(
                                    model, data, lambda s: fitting.initial(data, spec, s), c["fit_seeds"], folder,
                                    paths=c["pathfinder_paths"], draws=c["pathfinder_draws"],
                                    minimum_ess=c["minimum_importance_ess"], maximum_k=c["maximum_pareto_k"],
                                    keep_csv=c.get("keep_stan_csv", False))
                            result = dict(candidate=cand, variant=variant, **result)
                            if result["status"] == "ok":
                                result.update(summaries(sv, weights, obs, cand, c, scenario))
                            result.update(scenario=scenario, pool=pool, replicate=rep, source=source,
                                          selected_blocks=selected)
                            save_json(result_path, result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["plan", "compile", "simulate", "fit", "summarize", "calibrate"])
    p.add_argument("--config", type=Path, default=HERE / "config.json")
    p.add_argument("--out", type=Path, default=HERE / "runs" / "default")
    p.add_argument("--scenarios", nargs="+", choices=SCENARIOS, default=SCENARIOS)
    p.add_argument("--sources", nargs="+", choices=SOURCES, default=SOURCES)
    p.add_argument("--variants", nargs="+", choices=list(VARIANTS), default=list(VARIANTS))
    p.add_argument("--hapibd-command", type=lambda s: s.split(), default=None,
                   help="argv prefix, e.g. 'java -jar /path/hap-ibd.jar' (default: $HAPIBD_COMMAND or genetics_env)")
    p.add_argument("--block-limit", type=int, help="Debug only: prepare first N blocks")
    p.add_argument("--replicate-limit", type=int, help="Debug/pilot: fit first N subsamples")
    p.add_argument("--only-graphs", type=int, nargs="+", help="Debug only: partial rankings remain explicitly incomplete")
    p.add_argument("--candidate-slice", help="k/n: this worker fits candidates with index % n == k")
    p.add_argument("--block-slice", help="k/n: this worker simulates blocks with index % n == k")
    p.add_argument("--retry-failed", action="store_true")
    p.add_argument("--warmup", type=int, default=1000)
    p.add_argument("--samples", type=int, default=1000)
    p.add_argument("--parallel-chains", type=int, default=4)
    a = p.parse_args(); a.out = a.out.resolve()
    c = load_config(a.config)
    m = open_run(c, a.out)
    if a.stage == "plan":
        xs = m["candidate_orders"]
        nfits = len(xs)*len(a.variants)*len(a.sources)*len(a.scenarios)*c["pools"]*c["replicates"]
        print(f"{len(set(x['graph'] for x in xs))} graph shapes, {len(xs)} valid event orders")
        print(f"{c['blocks']} x {c['block_cm']} cM per pool; {c['genome_cm']} cM per subsample; no held-out blocks")
        print(f"{c['replicates']} overlapping subsamples/pool; {nfits} candidate fits (each with {len(c['fit_seeds'])} Pathfinder starts)")
        print(f"Correct backbone: {[x['id'] for x in xs if x['correct']]}")
        print(f"Manifest: {a.out / 'manifest.json'}")
    elif a.stage == "compile":
        # Build everything workers share BEFORE they start, so eight processes
        # never race on the same executable or the MRCA scanner's .so file.
        ibd.mrca_library()
        for variant in a.variants:
            model = models.normalized_model(VARIANTS[variant], a.out / "models")
            print(f"{variant}: {model.exe_file}")
    elif a.stage == "simulate":
        simulate(c, a)
    elif a.stage in ("fit", "calibrate"):
        fit_stage(c, a, calibration=a.stage == "calibrate")
    else:
        from visualize import build_summary
        build_summary(c, a.out)


if __name__ == "__main__":
    main()
