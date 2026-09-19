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
                   UPPER, pair_counts as data_pair_counts, load_config, edges, candidates, selections, simulation_demography,
                   quiet, to_msprime_demography, true_ibd, hapibd, snp_summaries,
                   aggregate, stan_data, save_json, digest, event_parameters)


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
    mapping = event_parameters(candidate, c)
    if not mapping:
        return None
    times = [x["truth"] for x in mapping if x["variable"] == "cumulative_times"]
    return dict(times=times, correct_order=bool(np.all(np.diff(times) > 0)))


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
    result["semantic_parameters"] = {}
    for item in event_parameters(candidate, c):
        values = sv[item["variable"]][:, item["index"]]
        if item["fold"]:
            values = np.maximum(values, 1-values)
        stat = parameter_summary(values, weights)
        stat.update(truth=item["truth"], variable=item["variable"], index=item["index"], folded=item["fold"])
        result["semantic_parameters"][item["name"]] = stat
    expected = np.maximum(pred["ibd_number"] * obs["cm"] * (np.asarray(data_pair_counts(c))), 1e-12)
    counts = obs["ibd_count"]
    from scipy.special import xlogy
    deviance = 2 * (xlogy(counts, counts / expected) - counts + expected)
    result["residuals"] = {
        "ibd_fraction": residual.tolist(),
        "ibd_count": (counts-expected).tolist(),
        "ibd_pearson": ((counts-expected)/np.sqrt(expected)).tolist(),
        "ibd_deviance": (np.sign(counts-expected)*np.sqrt(np.maximum(deviance, 0))).tolist(),
        "snp_raw": (obs["w_hat"]-pred["W_centered"]).tolist(),
        "snp_standardized": ((obs["w_hat"]-pred["W_centered"])/obs["w_se"]).tolist(),
    }
    result["observations"] = {k: np.asarray(v).tolist() for k,v in obs.items()}
    result["prediction_intervals"] = {name: parameter_summary(sv[name], weights)
                                       for name in ("ibd_fraction", "ibd_number", "W_centered")}
    result["likelihood_summaries"] = {name: parameter_summary(sv[name], weights)
        for name in ("lp_ibd", "lp_snp", "chi2_ibd", "chi2_snp") if name in sv}
    result["ne_log_distortion"] = {k: np.log(np.asarray(result[k]["mean"])/c["haploid_ne"]).tolist()
                                    for k in ("Ne", "Ne_ibd", "Ne_snp") if k in result}
    if "Ne_ibd" in result:
        result["log_ne_ibd_over_snp"] = parameter_summary(np.log(sv["Ne_ibd"]/sv["Ne_snp"]), weights)
    return result


def fit_one(model, data, obs, c, variant, candidate, folder):
    import arviz as az
    retained = defaultdict(list); logs = []; runs = []; total_proposal_draws = 0
    needed = ("cumulative_times", "admixture_fractions", "Ne", "Ne_ibd", "Ne_snp",
              "ibd_fraction", "ibd_number", "W_centered", "times",
              "lp_ibd", "lp_snp", "chi2_ibd", "chi2_snp")
    t0 = time.monotonic()
    for seed in c["fit_seeds"]:
        run = {"seed": seed}
        try:
            fit = model.pathfinder(data=data, inits=initial(data, variant, seed), seed=seed,
                  num_paths=c["pathfinder_paths"], draws=c["pathfinder_draws"],
                  num_single_draws=c["pathfinder_draws"], psis_resample=False, calculate_lp=True,
                  output_dir=str(folder / f"seed_{seed}"), sig_figs=12, show_console=False)
            names = list(fit.column_names)
            arr = np.asarray(fit.draws()).reshape(-1, len(names))
            lw = arr[:, names.index("lp__")] - arr[:, names.index("lp_approx__")]
            finite = np.isfinite(lw)
            if np.sum(finite) < 20:
                raise ValueError("Fewer than 20 finite importance weights")
            total_proposal_draws += len(lw)
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
                if name in needed or name.startswith(("mu_log", "sigma_log", "tau", "Ne_raw")):
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
                logz=float(logsumexp(lw)-np.log(total_proposal_draws)), elbo=float(lw.mean()),
                elbo_finite_draws_only=any(r.get("discarded", 0) for r in runs))
    sv = {name: np.concatenate(values) for name, values in retained.items()}
    # Keep draws of ALL sampled parameters plus sizes/times and component scores.
    # Prediction draws can be regenerated with Stan from these parameters.
    archive = {k: v for k,v in sv.items() if k not in ("ibd_fraction", "ibd_number", "W_centered")}
    np.savez_compressed(folder / "posterior_draws.npz", log_weight=lw, weight=weights, **archive)
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
                            if scenario == "no_loop":
                                for name, stat in result.get("semantic_parameters", {}).items():
                                    if name.startswith("loop_"):
                                        stat["truth"] = None
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
    from visualize import build_summary
    build_summary(c, a.out)


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
