#!/usr/bin/env python3
"""Resumable stages: plan, compile, simulate, fit (Pathfinder), nuts (NUTS + bridge), summarize.

`simulate --block-slice k/n`, `fit --task-slice k/n` and `nuts --task-slice k/n` let
several workers share one run directory; each result lives in its own folder and is
skipped once written.  `nuts` starts every chain from the highest-weight draws of the
Pathfinder fit of the same subsample, source, model and graph, so run `fit` first.

A block cache is complete once it holds every configured IBD source: `simulate`
adds only the sources a cache lacks (hap-IBD to a true-IBD cache, say), reading the
tree sequence back rather than simulating it again.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from study import (HERE, ROOT, POPS, GRAPHS, MODEL_NAMES, load_config, edges, n_haploid, parameters,
                   selections, fit_tasks, nuts_tasks, msprime_demography, true_ibd, hapibd, snp_summaries, aggregate,
                   stan_data, save_json)
from methods import evidence, fitting, ibd, models, parallel, simulate as sim
from methods.blocks import load_blocks
from methods.diagnostics import fit_summary
from methods.fitting import parameter_summary
from methods.utils import file_hashes, open_run


def manifest(c):
    files = [HERE / "study.py", HERE / "run_study.py", *sorted((ROOT / "methods").glob("*.py")),
             ROOT / "methods" / "mrca_scan.cpp", *[spec.path for spec in MODEL_NAMES.values()]]
    return {"config": c, "code_sha256": file_hashes(files, ROOT),
            "graphs": {g: build().ordered_events for g, build in GRAPHS.items()},
            "selected_blocks": {str(cm): selections(c, 0, cm) for cm in c["genome_cm"]}}


def pool(out):
    return out / "pool_000"


def result_folder(out, stage, cm, rep, source, model, graph):
    return pool(out) / stage / f"cm_{int(cm):04d}" / f"rep_{rep:03d}" / source / model / graph


def hapibd_provenance(command):
    prov = {"hapibd_command": command, "hapibd_version_output": ibd.hapibd_version(command)[:5000]}
    for item in command:
        if item.endswith(".jar") and Path(item).is_file():
            prov["jar_sha256"] = hashlib.sha256(Path(item).read_bytes()).hexdigest()
    return prov


def simulate(c, a):
    """Prepare the block caches.  Simulate every block, or -- when the config names
    `trees_from` -- re-read the tree sequences an earlier run already simulated; then
    extract SNP sums and every configured IBD source a cache still lacks.  hap-IBD's
    input VCF is deleted after the call (it is regenerated from the trees on demand);
    its segment files and log stay beside the cache."""
    import msprime
    import tskit
    folder = pool(a.out); folder.mkdir(parents=True, exist_ok=True)
    source = (ROOT / c["trees_from"]).resolve() if "trees_from" in c else None
    provenance = {"msprime": msprime.__version__, "tskit": tskit.__version__, "numpy": np.__version__,
                  "trees_from": str(source) if source else None}
    command = ibd.hapibd_command() if "hapibd" in c["sources"] else None
    if command:
        hap = hapibd_provenance(command)
        saved = folder / "provenance_hapibd.json"
        if saved.exists() and json.loads(saved.read_text()) != hap:
            raise ValueError("hap-IBD command or jar changed; use a new run directory.")
        save_json(saved, hap)
    save_json(folder / "provenance.json", provenance)
    msdem = None if source else msprime_demography(c)
    for block in range(c["blocks"]):
        if not parallel.in_slice(block, a.block_slice):
            continue
        cache = folder / f"block_{block:03d}.npz"
        values = dict(np.load(cache)) if cache.exists() else {}
        missing = [s for s in c["sources"] if s + "_count" not in values]
        if not missing:
            continue
        t0 = time.monotonic()
        trees = (source or folder) / f"block_{block:03d}.trees"
        if trees.exists():
            ts = tskit.load(trees)
            seed = int(np.load(source / f"block_{block:03d}.npz")["seed"]) if source else sim.block_seed(c["seed"], 0, block)
        else:
            seed = sim.block_seed(c["seed"], 0, block)
            ts = sim.simulate_block(msdem, {p: c["diploid_samples"] for p in POPS}, c["block_cm"],
                                    c["recombination_rate"], c["mutation_rate"], seed, model=c["ancestry_model"])
            ts.dump(trees)
        if "snp_numer" not in values:
            values.update(snp_summaries(ts, c))
        for src in missing:
            if src == "true":
                count, length = true_ibd(ts, c)
            else:
                work = folder / f"hap_{block:03d}"
                count, length = hapibd(ts, c, work, command)
                (work / "input.vcf").unlink()
            values[src + "_count"], values[src + "_length"] = count, length
        values["seed"] = np.asarray(seed)
        tmp = cache.with_suffix(".tmp.npz"); np.savez_compressed(tmp, **values); tmp.replace(cache)
        print(f"block {block}: {ts.num_trees} trees, {int(values['snp_sites'].sum())} SNPs, segments "
              + ", ".join(f"{s} {int(values[s + '_count'].sum())}" for s in c["sources"])
              + f" (added {', '.join(missing)}); {time.monotonic() - t0:.0f}s", flush=True)


def recovery(sv, w, c):
    return {p["name"]: dict(parameter_summary(sv[p["variable"]][:, p["index"]], w), truth=p["truth"])
            for p in parameters(c)}


def fit_stage(c, a):
    work = [t for i, t in enumerate(fit_tasks(c, a.sources, a.models, a.graphs, a.replicate_limit))
            if parallel.in_slice(i, a.task_slice)]
    compiled = {m: models.normalized_model(MODEL_NAMES[m], a.out / "models") for m in a.models}
    blocks = None
    for cm, rep, source, model, graph in work:
        out = result_folder(a.out, "pathfinder", cm, rep, source, model, graph)
        if (out / "result.json").exists():
            continue
        blocks = blocks or load_blocks(pool(a.out), c["blocks"], a.sources)
        selected = selections(c, 0, cm)[rep]
        obs = aggregate(blocks, selected, c, source)
        spec = MODEL_NAMES[model]; data = stan_data(graph, obs, c, model)
        out.mkdir(parents=True, exist_ok=True)
        print(f"cm={cm} rep={rep} {source} {model} {graph}", flush=True)
        result, sv, w = fitting.pathfinder(compiled[model], data, lambda s: fitting.initial(data, spec, s),
                                           c["fit_seeds"], out, paths=c["pathfinder_paths"],
                                           draws=c["pathfinder_draws"], minimum_ess=c["minimum_importance_ess"],
                                           maximum_k=c["maximum_pareto_k"], keep_csv=c["keep_stan_csv"])
        result = dict(model=model, graph=graph, **result)
        if result["status"] == "ok":
            result.update(fit_summary(sv, w, obs, n_haploid(c), edges(c)))
            if graph == "T_true":
                result["semantic_parameters"] = recovery(sv, w, c)
        result.update(genome_cm=cm, replicate=rep, source=source, selected_blocks=selected)
        save_json(out / "result.json", result)


def nuts_inits(out, cm, rep, source, model, graph, k):
    z = np.load(result_folder(out, "pathfinder", cm, rep, source, model, graph) / "posterior_draws.npz")
    top = np.argsort(-z["log_weight"])[:k]
    names = [n for n in z.files if n.startswith(("times", "admixture_fractions", "mu_log", "sigma_log", "tau", "Ne_raw"))]
    return [{n: z[n][i].tolist() for n in names} for i in top]


def predictive(sv, obs, c, seed=0):
    """Posterior-predictive IBD counts per bin and pair: Poisson(lambda) for every draw."""
    lam = np.asarray(sv["ibd_number"]) * obs["cm"] * ibd.pair_counts(n_haploid(c))
    counts = np.random.default_rng(seed).poisson(np.maximum(lam, 0))
    q = np.percentile(counts, [5, 50, 95], axis=0)
    return {"q05": q[0].tolist(), "q50": q[1].tolist(), "q95": q[2].tolist(), "lambda_mean": lam.mean(0).tolist()}


def nuts_stage(c, a):
    work = [t for i, t in enumerate(nuts_tasks(c, a.graphs)) if parallel.in_slice(i, a.task_slice)]
    n = c["nuts"]
    compiled = models.normalized_model(MODEL_NAMES["mixed"], a.out / "models")
    blocks = None
    for cm, rep, source, model, graph in work:
        out = result_folder(a.out, "nuts", cm, rep, source, model, graph)
        if (out / "result.json").exists():
            continue
        blocks = blocks or load_blocks(pool(a.out), c["blocks"], [source])
        obs = aggregate(blocks, selections(c, 0, cm)[rep], c, source)
        spec = MODEL_NAMES[model]; data = stan_data(graph, obs, c, model)
        out.mkdir(parents=True, exist_ok=True)
        print(f"cm={cm} rep={rep} {source} {model} {graph} NUTS", flush=True)
        t0 = time.monotonic()
        result, sv, w = fitting.nuts(compiled, data, nuts_inits(a.out, cm, rep, source, model, graph, n["chains"]), c["seed"],
                                     out, chains=n["chains"], warmup=n["warmup"], samples=n["samples"],
                                     parallel_chains=a.parallel_chains, adapt_delta=n["adapt_delta"])
        result = dict(model=model, graph=graph, **result)
        if result["status"] == "ok":
            result["evidence"] = evidence.log_evidence(compiled, data, spec, sorted((out / "chains").glob("*.csv")),
                                                       out / "bridge", seed=c["seed"])
            result.update(fit_summary(sv, w, obs, n_haploid(c), edges(c)))
            result["posterior_predictive"] = predictive(sv, obs, c)
            if graph == "T_true":
                result["semantic_parameters"] = recovery(sv, w, c)
        result.update(genome_cm=cm, replicate=rep, source=source, seconds_total=time.monotonic() - t0)
        save_json(out / "result.json", result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["plan", "compile", "simulate", "fit", "nuts", "summarize"])
    p.add_argument("--config", type=Path, default=HERE / "config.json")
    p.add_argument("--out", type=Path, default=HERE / "runs" / "default")
    p.add_argument("--sources", nargs="+", help="IBD sources to fit (default: the config's)")
    p.add_argument("--models", nargs="+", choices=list(MODEL_NAMES), default=list(MODEL_NAMES))
    p.add_argument("--graphs", nargs="+", choices=list(GRAPHS), default=list(GRAPHS))
    p.add_argument("--replicate-limit", type=int)
    p.add_argument("--task-slice"); p.add_argument("--block-slice")
    p.add_argument("--parallel-chains", type=int, default=4)
    a = p.parse_args(); a.out = a.out.resolve()
    c = load_config(a.config)
    a.sources = a.sources or c["sources"]
    if not set(a.sources) <= set(c["sources"]):
        raise ValueError(f"--sources must be among the config's {c['sources']}")
    if a.stage == "summarize":
        saved = json.loads((a.out / "manifest.json").read_text())
        if saved["config"] != c:
            raise ValueError("Configuration differs from the run being summarised.")
        from visualize import build_summary
        build_summary(c, a.out)
        return
    open_run(a.out, manifest(c))
    if a.stage == "plan":
        print(f"{len(fit_tasks(c, a.sources, a.models, a.graphs))} Pathfinder fits ({', '.join(a.sources)} x "
              f"{', '.join(a.models)} x {', '.join(a.graphs)}); {len(nuts_tasks(c, a.graphs))} NUTS fits "
              f"({c['nuts_source']}, mixed x {', '.join(a.graphs)} x {c['nuts_replicates']} subsamples)")
        print(f"{c['blocks']} x {c['block_cm']} cM blocks; {len(c['bin_edges_cm']) - 1} bins; lengths {c['genome_cm']}")
    elif a.stage == "compile":
        ibd.mrca_library()
        for m in MODEL_NAMES:
            print(m, models.normalized_model(MODEL_NAMES[m], a.out / "models").exe_file)
    elif a.stage == "simulate":
        simulate(c, a)
    elif a.stage == "fit":
        fit_stage(c, a)
    else:
        nuts_stage(c, a)


if __name__ == "__main__":
    main()
