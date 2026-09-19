#!/usr/bin/env python3
"""Resumable stages: plan, simulate, fit, summarize, and calibrate (MCMC)."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import sys
import time

import numpy as np
from scipy.special import logsumexp

from study import (HERE, PIPELINE, REAL_DATA, POPS, SCENARIOS, SOURCES, VARIANTS,
                   UPPER, load_config, edges, candidates, selections, simulation_demography,
                   quiet, to_msprime_demography, true_ibd, hapibd, snp_summaries,
                   aggregate, stan_data, save_json, digest)


def manifest(c):
    files = [HERE / "study.py", HERE / "run_study.py", HERE / "mrca_scan.cpp",
             PIPELINE / "demography.py", PIPELINE / "simulation_methods.py",
             REAL_DATA / "3pop" / "enumerate_3pop.py"]
    files += [REAL_DATA / f for f in VARIANTS.values()]
    return {"config": c, "code_sha256": {str(p.relative_to(PIPELINE)):
            hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            "candidate_orders": candidates(),
            "selected_blocks": {str(p): selections(c, p) for p in range(c["pools"])}}


def open_run(c, path):
    path.mkdir(parents=True, exist_ok=True)
    new = manifest(c)
    saved = path / "manifest.json"
    if saved.exists():
        old = json.loads(saved.read_text())
        if digest(old) != digest(new):
            raise ValueError("Configuration or source code changed. Use a new --out directory; cached results must not be mixed.")
    else:
        save_json(saved, new)
    return new


def pool_folder(out, scenario, pool):
    return out / scenario / f"pool_{pool:03d}"


def simulate(c, a):
    import msprime
    import tskit
    command = shlex.split(a.hapibd_command) if a.hapibd_command else ["hap-ibd"]
    if "hapibd" in a.sources and not shutil.which(command[0]):
        raise RuntimeError("hap-IBD launcher missing. Set --hapibd-command 'java -jar /path/hap-ibd.jar'.")
    provenance = {"msprime": msprime.__version__, "tskit": tskit.__version__,
                  "numpy": np.__version__, "hapibd_command": command if "hapibd" in a.sources else None}
    if "hapibd" in a.sources:
        # Fail before ancestry simulation if Java/JAR cannot be launched.
        check = subprocess_run_version(command)
        provenance["hapibd_version_output"] = check[:5000]
        for item in command:
            if item.endswith(".jar") and Path(item).is_file():
                provenance["jar_sha256"] = hashlib.sha256(Path(item).read_bytes()).hexdigest()
    for scenario in a.scenarios:
        dem = quiet(simulation_demography, c, scenario)
        msdem = quiet(to_msprime_demography, dem)
        msdem.validate()
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
                cache = folder / f"block_{block:03d}.npz"
                values = dict(np.load(cache)) if cache.exists() else {}
                if all(source+"_count" in values for source in a.sources):
                    continue
                seed = int(np.random.SeedSequence([c["seed"], SCENARIOS.index(scenario), pool, block]).generate_state(1)[0])
                seed = seed or 1
                tree_file = folder / f"block_{block:03d}.trees"
                t0 = time.monotonic()
                print(f"{scenario} pool {pool} block {block+1}/{c['blocks']}: ancestry/SNPs", flush=True)
                if tree_file.exists():
                    ts = tskit.load(tree_file)
                else:
                    ancestry = msprime.sim_ancestry(samples={p: c["diploid_samples"] for p in POPS},
                        ploidy=2, demography=msdem,
                        sequence_length=round(c["block_cm"]/(100*c["recombination_rate"])),
                        recombination_rate=c["recombination_rate"], random_seed=seed,
                        model=c["ancestry_model"], discrete_genome=True)
                    ts = msprime.sim_mutations(ancestry, rate=c["mutation_rate"], random_seed=seed,
                                              discrete_genome=True)
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


def subprocess_run_version(command):
    import subprocess
    p = subprocess.run(command, capture_output=True, text=True, timeout=30)
    text = p.stdout + p.stderr
    # hap-IBD prints usage without arguments; allow its usage exit status.
    if "hap-ibd" not in text.lower() or "Unable to locate a Java Runtime" in text:
        raise RuntimeError(f"Could not launch hap-IBD: {text[-2000:]}")
    return text


def normalized_model(variant, out):
    """Keep original model behavior; restore all density constants for logZ.

    Stan drops constants in ~ statements. Replacing them with explicit lpdf
    also normalizes the truncated exponential and half-Normal priors. This
    changes only additive constants, never posterior distributions.
    """
    from cmdstanpy import CmdStanModel
    source = (REAL_DATA / VARIANTS[variant]).read_text()
    pattern = r"(?m)^(\s*)([\w]+(?:\[[^\]\n]+\])?)\s*~\s*(normal|std_normal|exponential|beta)\(([^;\n]*)\);"
    def replace(m):
        indent, lhs, distribution, args = m.groups()
        return f"{indent}target += {distribution}_lpdf({lhs}" + (f" | {args}" if args.strip() else "") + ");"
    source, n = re.subn(pattern, replace, source)
    expected = 11 if variant.endswith("separate") else 7
    if n != expected:
        raise ValueError(f"Prior/likelihood source changed: expected {expected} sampling statements, found {n}")
    ntraj = 2 if variant.endswith("separate") else 1
    source = source.replace("\nmodel {", "\nmodel {\n    // Normalizers for times >= 1 and two half-Normals per trajectory.\n"
                            f"    target += n_events * 0.01 + {2*ntraj} * log(2.0);", 1)
    folder = out / "models"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / VARIANTS[variant]
    if not path.exists() or path.read_text() != source:
        path.write_text(source)
    return CmdStanModel(stan_file=str(path))


def initial(data, variant, seed):
    # Dispersed, truth-independent starts; same initialization schedule for every graph.
    rng = np.random.default_rng(seed)
    scale = (10., 50., 150., 400.)[seed % 4]
    init = {"times": np.maximum(1.01, scale*np.exp(rng.normal(0, 0.4, data["n_events"]))),
            "admixture_fractions": rng.uniform(0.1, 0.9, data["n_admixture"])}
    suffixes = ("_ibd", "_snp") if variant.endswith("separate") else ("",)
    level = np.log(15000.) + rng.normal(0, 0.4)
    for suffix in suffixes:
        init.update({"mu_log"+suffix: level, "sigma_log"+suffix: 0.2, "tau"+suffix: 0.2,
                     "Ne_raw"+suffix: rng.normal(0, 0.2, data["n_nodes"])})
    return init


def parameter_summary(draws, weights):
    draws = np.asarray(draws)
    shape = draws.shape[1:]
    flat = draws.reshape(len(weights), -1)
    mean = weights @ flat
    bounds = []
    for j in range(flat.shape[1]):
        order = np.argsort(flat[:, j])
        cdf = np.cumsum(weights[order])
        bounds.append(np.interp([0.025, 0.975], cdf, flat[order, j]))
    bounds = np.asarray(bounds).reshape(-1, 2)
    return {"mean": mean.reshape(shape).tolist(),
            "lower": bounds[:, 0].reshape(shape).tolist(),
            "upper": bounds[:, 1].reshape(shape).tolist()}


def truth_for(candidate, c):
    if not candidate["correct"]:
        return None
    # Match parameters by ancestry clade, never by event number or source label.
    clades = {name: {name} for name in POPS}
    times = []; a_source = None
    for ev in candidate["events"]:
        if ev["type"] == "ADMIXTURE":
            for parent in ev["parents"]:
                clades[parent] = {parent}
            times.append(c["times"]["b_split"])
        else:
            members = set.union(*(clades[x] for x in ev["children"]))
            clades[ev["parent"]] = members
            which = "root" if {"a", "c"} <= members else "left" if "a" in members else "right"
            times.append(c["times"][which])
            if which == "left":
                a_source = next(x for x in members if x.startswith("b."))
    adm = next(ev for ev in candidate["events"] if ev["type"] == "ADMIXTURE")
    fraction = c["fractions"]["b"] if adm["parents"][0] == a_source else 1-c["fractions"]["b"]
    return {"times": times, "fraction": fraction, "correct_order": bool(np.all(np.diff(times) > 0))}


def summaries(sv, weights, obs, candidate, c):
    result = {}
    for name in ("cumulative_times", "admixture_fractions", "Ne", "Ne_ibd", "Ne_snp"):
        if name in sv:
            result[name] = parameter_summary(sv[name], weights)
    pred = {name: np.tensordot(weights, sv[name], axes=1)
            for name in ("ibd_fraction", "ibd_number", "W_centered")}
    result["predictions"] = {name: value.tolist() for name, value in pred.items()}
    residual = obs["ibd_hat"] - pred["ibd_fraction"]
    result["fit_checks"] = {
        "ibd_rmse_by_pair": np.sqrt(np.mean(residual[:, UPPER[0], UPPER[1]]**2, axis=0)).tolist(),
        "ibd_residual_by_bin_pair": residual[:, UPPER[0], UPPER[1]].tolist(),
        "snp_rmse": float(np.sqrt(np.mean((obs["w_hat"][UPPER]-pred["W_centered"][UPPER])**2))),
        "snp_standardized_residual": ((obs["w_hat"]-pred["W_centered"])/obs["w_se"]).tolist()}
    truth = truth_for(candidate, c)
    if truth:
        result["truth"] = truth
        if truth["correct_order"]:
            ts = result["cumulative_times"]; fs = result["admixture_fractions"]
            result["recovery"] = {
                "time_relative_error": ((np.asarray(ts["mean"])-truth["times"])/truth["times"]).tolist(),
                "time_covered": ((np.asarray(ts["lower"]) <= truth["times"]) &
                                 (np.asarray(ts["upper"]) >= truth["times"])).tolist(),
                "time_width": (np.asarray(ts["upper"])-ts["lower"]).tolist(),
                "fraction_error": fs["mean"][0]-truth["fraction"],
                "fraction_covered": fs["lower"][0] <= truth["fraction"] <= fs["upper"][0],
                "fraction_width": fs["upper"][0]-fs["lower"][0]}
    result["ne_log_distortion"] = {k: np.log(np.asarray(result[k]["mean"])/c["haploid_ne"]).tolist()
                                    for k in ("Ne", "Ne_ibd", "Ne_snp") if k in result}
    if "Ne_ibd" in result:
        result["log_ne_ibd_over_snp"] = parameter_summary(np.log(sv["Ne_ibd"]/sv["Ne_snp"]), weights)
    return result


def fit_one(model, data, obs, c, variant, candidate, folder):
    import arviz as az
    retained = defaultdict(list); logs = []; runs = []
    needed = ("cumulative_times", "admixture_fractions", "Ne", "Ne_ibd", "Ne_snp",
              "ibd_fraction", "ibd_number", "W_centered")
    t0 = time.monotonic()
    for seed in c["fit_seeds"]:
        run = {"seed": seed}
        try:
            fit = model.pathfinder(data=data, inits=initial(data, variant, seed), seed=seed,
                  num_paths=c["pathfinder_paths"], draws=c["pathfinder_draws"],
                  num_single_draws=c["pathfinder_draws"], psis_resample=False, calculate_lp=True,
                  output_dir=str(folder / f"seed_{seed}"), show_console=False)
            names = list(fit.column_names)
            arr = np.asarray(fit.draws()).reshape(-1, len(names))
            lw = arr[:, names.index("lp__")] - arr[:, names.index("lp_approx__")]
            finite = np.isfinite(lw)
            if np.sum(finite) < 20:
                raise ValueError("Fewer than 20 finite importance weights")
            lw = lw[finite]
            _, k = az.psislw(lw)
            rawweights = np.exp(lw-logsumexp(lw))
            run.update(logz=float(logsumexp(lw)-np.log(len(lw))), elbo=float(lw.mean()),
                       ess=float(1/np.sum(rawweights**2)), pareto_k=float(k), draws=len(lw),
                       discarded=int(np.sum(~finite)))
            # Infinity is an explicit diagnostic failure, not invalid JSON.
            if not np.isfinite(run["pareto_k"]):
                run["pareto_k"] = None
            for name, values in fit.stan_variables().items():
                if name in needed:
                    retained[name].append(np.asarray(values)[finite])
            logs.append(lw)
        except Exception as ex:
            run["error"] = str(ex)
        runs.append(run)
    base = {"candidate": candidate, "variant": variant, "runs": runs,
            "seconds": time.monotonic()-t0, "method": "pathfinder_importance_sampling"}
    if not logs:
        return dict(base, status="failed", reliable=False)
    lw = np.concatenate(logs)
    np.savez_compressed(folder / "importance_weights.npz",
                        **{f"start_{i}": values for i, values in enumerate(logs)})
    weights = np.exp(lw-logsumexp(lw))
    _, k = az.psislw(lw)
    ess = float(1/np.sum(weights**2))
    reliable = (len(logs) == len(c["fit_seeds"]) and ess >= c["minimum_importance_ess"]
                and np.isfinite(k) and float(k) <= c["maximum_pareto_k"]
                and all(r.get("discarded", 1) == 0 for r in runs))
    base.update(status="ok", reliable=bool(reliable), ess=ess,
                pareto_k=float(k) if np.isfinite(k) else None,
                logz=float(logsumexp(lw)-np.log(len(lw))), elbo=float(lw.mean()))
    sv = {name: np.concatenate(values) for name, values in retained.items()}
    base.update(summaries(sv, weights, obs, candidate, c))
    return base


def load_blocks(folder, c, sources):
    blocks = []
    for k in range(c["blocks"]):
        path = folder / f"block_{k:03d}.npz"
        with np.load(path) as z:
            block = dict(z)
        if any(s+"_count" not in block for s in sources):
            raise ValueError(f"Missing IBD source in {path}; finish simulate first")
        blocks.append(block)
    return blocks


def fit_stage(c, a, calibration=False):
    chosen_candidates = candidates()
    if calibration:
        chosen_candidates = [x for x in chosen_candidates if x["correct"] and truth_for(x, c)["correct_order"]]
    if a.only_graphs:
        chosen_candidates = [x for x in chosen_candidates if x["graph"] in a.only_graphs]
    if not chosen_candidates:
        raise ValueError("No requested candidate graphs")
    for variant in a.variants:
        model = normalized_model(variant, a.out)
        for scenario in a.scenarios:
            for pool in range(c["pools"]):
                pf = pool_folder(a.out, scenario, pool)
                blocks = load_blocks(pf, c, a.sources)
                for rep, selected in enumerate(selections(c, pool)):
                    if a.replicate_limit is not None and rep >= a.replicate_limit:
                        break
                    for source in a.sources:
                        obs = aggregate(blocks, selected, source, c)
                        for cand in chosen_candidates:
                            folder = pf / ("mcmc" if calibration else "fits") / f"rep_{rep:03d}" / source / variant / cand["id"]
                            result_path = folder / "result.json"
                            if result_path.exists() and not a.retry_failed:
                                continue
                            if result_path.exists() and json.loads(result_path.read_text()).get("reliable"):
                                continue
                            folder.mkdir(parents=True, exist_ok=True)
                            data = stan_data(cand, obs, c)
                            print(f"{scenario} pool={pool} rep={rep} {source} {variant} {cand['id']}", flush=True)
                            if calibration:
                                result = mcmc_one(model, data, obs, c, variant, cand, folder, a)
                            else:
                                result = fit_one(model, data, obs, c, variant, cand, folder)
                            result.update(scenario=scenario, pool=pool, replicate=rep, source=source,
                                          selected_blocks=selected)
                            save_json(result_path, result)
                            # Keep compact scores/weights and logs. Large generated-quantity
                            # CSVs across 43,200 fits otherwise require terabytes of storage.
                            # Remove only our own Pathfinder CSVs, after the result is saved.
                            if not calibration and result["status"] == "ok" and not c.get("keep_stan_csv", False):
                                for path in folder.glob("seed_*/*.csv"):
                                    path.unlink()


def mcmc_one(model, data, obs, c, variant, candidate, folder, a):
    base = dict(candidate=candidate, variant=variant, method="mcmc",
                sampling_settings=dict(warmup=a.warmup, samples=a.samples, chains=4,
                                       parallel_chains=a.parallel_chains, adapt_delta=0.95))
    start = time.monotonic()
    try:
        fit = model.sample(data=data, inits=[initial(data, variant, seed) for seed in (1, 7, 13, 19)],
            seed=c["seed"], chains=4, parallel_chains=a.parallel_chains,
            iter_warmup=a.warmup, iter_sampling=a.samples, adapt_delta=0.95,
            output_dir=str(folder / "chains"), show_progress=False)
        diagnostics = fit.method_variables()
        table = fit.summary()
        # Constants/transformed matrices contain identically zero entries: omit
        # them from convergence summaries, but inspect all sampled parameters.
        keys = [k for k in table.index if re.match(r"^(times\[|admixture_fractions\[|mu_log|sigma_log|tau|Ne_raw)", k)]
        sub = table.loc[keys]
        rhat = float(sub["R_hat"].max()); ess = float(sub["ESS_bulk"].min())
        divergent = int(diagnostics["divergent__"].sum())
        depth_hits = int((diagnostics["treedepth__"] >= 10).sum())
        base.update(status="ok", reliable=bool(np.isfinite(rhat) and rhat <= 1.01 and ess >= 400
                                               and divergent == 0 and depth_hits == 0),
                    max_rhat=rhat if np.isfinite(rhat) else None,
                    min_ess_bulk=ess if np.isfinite(ess) else None,
                    divergences=divergent, max_treedepth_hits=depth_hits)
        sv = fit.stan_variables()
        weights = np.full(len(sv["cumulative_times"]), 1/len(sv["cumulative_times"]))
        base.update(summaries(sv, weights, obs, candidate, c))
    except Exception as ex:
        base.update(status="failed", reliable=False, error=str(ex))
    base["seconds"] = time.monotonic()-start
    return base


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        if rows:
            writer = csv.DictWriter(handle, list(rows[0]))
            writer.writeheader(); writer.writerows(rows)


def summarize(c, a):
    all_candidates = candidates()
    by_graph = defaultdict(list)
    for cand in all_candidates:
        by_graph[cand["graph"]].append(cand)
    records = []; rankings = []; recovery = []; fitrows = []
    for scenario in a.scenarios:
        for pool in range(c["pools"]):
            pf = pool_folder(a.out, scenario, pool)
            for source in a.sources:
                for variant in a.variants:
                    for rep in range(c["replicates"]):
                        folder = pf / "fits" / f"rep_{rep:03d}" / source / variant
                        fitted = {p.parent.name: json.loads(p.read_text()) for p in folder.glob("*/result.json")}
                        identity = dict(scenario=scenario, pool=pool, replicate=rep, source=source, variant=variant)
                        for cand_id, r in fitted.items():
                            fitrows.append(dict(identity, candidate=cand_id, status=r["status"], reliable=r["reliable"],
                                                seconds=r["seconds"], ess=r.get("ess"), pareto_k=r.get("pareto_k")))
                            if "recovery" in r:
                                recovery.append(dict(identity, method="pathfinder", reliable=r["reliable"], **r["recovery"]))
                        scores = []; complete = True; reliable = True
                        for graph, orders in by_graph.items():
                            fits = [fitted.get(x["id"]) for x in orders]
                            if any(r is None or r["status"] != "ok" for r in fits):
                                complete = False
                                continue
                            reliable &= all(r["reliable"] for r in fits)
                            # Equal prior over graph shapes; equal conditional prior over orders.
                            # Do not maximize across orders, which rewards shapes with more orders.
                            scores.append(dict(graph=graph, newick=orders[0]["newick"], admixed=orders[0]["admixed"],
                                               correct=orders[0]["correct"],
                                               logz=float(logsumexp([r["logz"] for r in fits])-np.log(len(fits))),
                                               elbo=float(logsumexp([r["elbo"] for r in fits])-np.log(len(fits)))))
                        row = dict(identity, complete=complete, reliable=bool(complete and reliable),
                                   winner=None, correct=None, correct_rank=None, score_gap=None,
                                   b_is_admixed=None, elbo_winner=None, elbo_correct=None)
                        if complete:
                            scores.sort(key=lambda x: x["logz"], reverse=True)
                            true = next(x for x in scores if x["correct"])
                            rival = max(x["logz"] for x in scores if not x["correct"])
                            ebest = max(scores, key=lambda x: x["elbo"])
                            row.update(winner=scores[0]["graph"], correct=scores[0]["correct"],
                                       correct_rank=1+sum(x["logz"] > true["logz"] for x in scores),
                                       score_gap=true["logz"]-rival, b_is_admixed=scores[0]["admixed"] == "b",
                                       elbo_winner=ebest["graph"], elbo_correct=ebest["correct"])
                            for rank, score in enumerate(scores, 1):
                                rankings.append(dict(identity, rank=rank, **score))
                        records.append(row)
                    for p in (pf / "mcmc").glob(f"rep_*/{source}/{variant}/*/result.json"):
                        r = json.loads(p.read_text())
                        if "recovery" in r:
                            recovery.append({k: r[k] for k in ("scenario", "pool", "replicate", "source", "variant")}
                                            | dict(method="mcmc", reliable=r["reliable"], **r["recovery"]))
    output = a.out / "summary"
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "topology_recovery.csv", records)
    write_csv(output / "topology_rankings.csv", rankings)
    write_csv(output / "fit_diagnostics.csv", fitrows)
    save_json(output / "parameter_recovery.json", recovery)
    groups = defaultdict(list)
    for row in records:
        groups[(row["scenario"], row["source"], row["variant"])].append(row)
    lines = ["# Hidden-loop topology robustness", "",
             "Recovery rates describe overlapping 30-of-50 block subsamples, not independent simulated genomes.",
             "One separate 2500 cM pool per scenario by default. No validation/test split.",
             "No binomial error bars or independent-replicate significance tests are used.", "",
             "Topology evidence is compared ONLY within each model variant and IBD source. Graph scores average",
             "over valid event orders with equal conditional prior weight. Pathfinder estimates are approximate.",
             "A replicate is reliable only when every candidate order passes the importance-weight checks.", "",
             "| Scenario | IBD | Model | Complete / requested | Reliable | Recovery (complete) | Recovery (reliable) |",
             "|---|---|---|---:|---:|---:|---:|"]
    def rate(rows):
        return f"{np.mean([r['correct'] for r in rows]):.1%}" if rows else "—"
    for key, rows in groups.items():
        complete = [r for r in rows if r["complete"]]
        reliable = [r for r in rows if r["reliable"]]
        lines.append("| " + " | ".join(key) + f" | {len(complete)} / {len(rows)} | {len(reliable)} | {rate(complete)} | {rate(reliable)} |")
    # Direct paired shared/separate comparison: missing or unreliable results never become failures of biology.
    paired = []
    lookup = {(r["scenario"], r["pool"], r["replicate"], r["source"], r["variant"]): r for r in records}
    for r in records:
        if not r["variant"].endswith("shared"):
            continue
        other = lookup.get((r["scenario"], r["pool"], r["replicate"], r["source"], r["variant"].replace("shared", "separate")))
        if other and r["complete"] and other["complete"]:
            paired.append({k: r[k] for k in ("scenario", "pool", "replicate", "source")} | dict(
                likelihood=r["variant"].split("_")[0], both_reliable=r["reliable"] and other["reliable"],
                shared_correct=r["correct"], separate_correct=other["correct"],
                recovery_difference=int(other["correct"])-int(r["correct"]),
                rank_improvement=r["correct_rank"]-other["correct_rank"]))
    write_csv(output / "paired_ne_comparison.csv", paired)
    parameter_rows = []
    pg = defaultdict(list)
    for r in recovery:
        if r["reliable"]:
            pg[tuple(r[k] for k in ("scenario", "source", "variant", "method"))].append(r)
    for key, rows in pg.items():
        for j, param in enumerate(("b_admixture_time", "left_merge_time", "right_merge_time", "root_time", "b_fraction")):
            err = np.asarray([r["time_relative_error"][j] if j < 4 else r["fraction_error"] for r in rows])
            cov = [r["time_covered"][j] if j < 4 else r["fraction_covered"] for r in rows]
            width = [r["time_width"][j] if j < 4 else r["fraction_width"] for r in rows]
            parameter_rows.append(dict(zip(("scenario", "source", "variant", "method"), key)) | dict(
                parameter=param, n=len(rows), bias=float(err.mean()), rmse=float(np.sqrt(np.mean(err**2))),
                coverage=float(np.mean(cov)), mean_interval_width=float(np.mean(width))))
    write_csv(output / "parameter_metrics.csv", parameter_rows)
    lines += ["", "Parameter metrics condition on fitting the correct backbone AND event order; they are not",
              "post-topology-selection coverage. Time errors are relative; fraction errors are absolute.",
              "Pathfinder interval coverage is approximate; use `calibrate` for MCMC coverage and diagnostics.",
              "With one pool, both coverage and bias summarize subsample behavior, not repeated independent simulations.",
              "", "Raw residual matrices, Ne distortions, separate-Ne contrasts, and per-start diagnostics are",
              "saved in each fit's result.json. The fitted-data residuals are descriptive, not held-out prediction."]
    (output / "report.md").write_text("\n".join(lines)+"\n")
    plot_results(records, output)
    print(f"Wrote {output / 'report.md'}", flush=True)


def plot_results(records, folder):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    variants = list(VARIANTS)
    for ax, scenario in zip(axes, SCENARIOS):
        for shift, source in zip((-0.18, 0.18), SOURCES):
            for j, variant in enumerate(variants):
                rows = [r for r in records if r["scenario"] == scenario and r["source"] == source
                        and r["variant"] == variant and r["reliable"]]
                if rows:
                    ax.bar(j+shift, np.mean([r["correct"] for r in rows]), width=0.34,
                           color="C0" if source == "true" else "C1", label=source if j == 0 else None)
                    ax.text(j+shift, 0.02, f"n={len(rows)}", ha="center", fontsize=7, rotation=90)
        ax.set_title(scenario); ax.set_xticks(range(4), [v.replace("_", "\n") for v in variants]); ax.set_ylim(0, 1.05)
        if ax.get_legend_handles_labels()[0]:
            ax.legend()
    axes[0].set_ylabel("Correct backbone frequency (reliable subsamples)")
    fig.suptitle("Overlapping subsamples; blank bars mean no reliable complete comparisons")
    fig.tight_layout(); fig.savefig(folder / "topology_recovery.png", dpi=160); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["plan", "simulate", "fit", "summarize", "calibrate"])
    p.add_argument("--config", type=Path, default=HERE / "config.json")
    p.add_argument("--out", type=Path, default=HERE / "runs" / "default")
    p.add_argument("--scenarios", nargs="+", choices=SCENARIOS, default=SCENARIOS)
    p.add_argument("--sources", nargs="+", choices=SOURCES, default=SOURCES)
    p.add_argument("--variants", nargs="+", choices=list(VARIANTS), default=list(VARIANTS))
    p.add_argument("--hapibd-command", default=os.environ.get("HAPIBD_COMMAND"), help="argv prefix, e.g. 'java -jar /path/hap-ibd.jar'")
    p.add_argument("--block-limit", type=int, help="Debug only: prepare first N blocks")
    p.add_argument("--replicate-limit", type=int, help="Debug/pilot: fit first N subsamples")
    p.add_argument("--only-graphs", type=int, nargs="+", help="Debug only: partial rankings remain explicitly incomplete")
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
    elif a.stage == "simulate":
        simulate(c, a)
    elif a.stage in ("fit", "calibrate"):
        fit_stage(c, a, calibration=a.stage == "calibrate")
    else:
        summarize(c, a)


if __name__ == "__main__":
    main()
