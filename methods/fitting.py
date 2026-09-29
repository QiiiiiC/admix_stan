"""Fitting engines: Pathfinder + importance sampling, and NUTS for calibration.

Pathfinder is fast and mode-seeking.  Every fit uses several dispersed,
truth-independent starts and pools their draws; the log importance weights
lw = lp__ - lp_approx__ give

  elbo  mean(lw)                            low-variance lower bound on log evidence
  elbo_max  best single path's mean(lw)     the tightest of the per-path lower bounds
  logz  logsumexp(lw) - log(#draws)         importance-sampled log evidence
  ess   1 / sum(w^2),  pareto_k from PSIS   whether logz can be trusted

In practice Pareto k > 1 is common on these models (the weights then have
infinite mean and logz is one draw's weight), so studies rank by ELBO and
report logz/ESS/k honestly.  Evidence is only comparable across graphs when the
model is normalised (`models.normalized_model`).
"""
from __future__ import annotations

from collections import defaultdict
import re
import time

import numpy as np
from scipy.special import logsumexp

KEEP = ("cumulative_times", "admixture_fractions", "times", "effective_N", "Ne", "Ne_ibd", "Ne_snp",
        "ibd_fraction", "ibd_number", "W_centered", "lp_ibd", "lp_snp", "chi2_ibd", "chi2_snp")
KEEP_PREFIX = ("mu_log", "sigma_log", "tau", "Ne_raw")
PREDICTIONS = ("ibd_fraction", "ibd_number", "W_centered")


def initial(data, spec, seed, ne_level=15000.0):
    """Dispersed, truth-independent start; the same schedule for every graph.

    Seeds cycle the time scale through 10/50/150/400 generations per increment
    so the starts cover recent-to-ancient solutions.
    """
    rng = np.random.default_rng(seed)
    scale = (10., 50., 150., 400.)[seed % 4]
    init = {"times": np.maximum(1.01, scale * np.exp(rng.normal(0, 0.4, data["n_events"]))),
            "admixture_fractions": rng.uniform(0.1, 0.9, data["n_admixture"])}
    level = np.log(ne_level) + rng.normal(0, 0.4)
    if spec.ne == "fixed":
        init["effective_N"] = float(np.exp(level))
    else:
        suffixes = ("_ibd", "_snp") if spec.trajectories == 2 else ("",)
        for suffix in suffixes:
            init.update({"mu_log" + suffix: level, "sigma_log" + suffix: 0.2,
                         "Ne_raw" + suffix: rng.normal(0, 0.2, data["n_nodes"])})
            if spec.ne == "smooth":
                init["tau" + suffix] = 0.2
    return init


def parameter_summary(draws, weights):
    """Weighted mean and 95% interval of every element of a draws array."""
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


def _kept(name):
    return name in KEEP or name.startswith(KEEP_PREFIX)


def pathfinder(model, data, inits, seeds, folder, paths=4, draws=1000,
               minimum_ess=100.0, maximum_k=0.7, keep_csv=False):
    """Multi-start Pathfinder with pooled importance weights.

    inits(seed) -> init dict.  Returns (result, draws, weights): `result` is
    JSON-ready (per-start and pooled elbo/logz/ess/k, `reliable`); `draws` holds
    the retained Stan variables stacked over starts; `weights` are normalised
    importance weights over those draws.  A start that fails is recorded with its
    error; the fit fails only if every start does.
    """
    import arviz as az
    retained = defaultdict(list); logs = []; runs = []; total = 0
    t0 = time.monotonic()
    for seed in seeds:
        run = {"seed": seed}
        try:
            fit = model.pathfinder(data=data, inits=inits(seed), seed=seed, num_paths=paths, draws=draws,
                                   num_single_draws=draws, psis_resample=False, calculate_lp=True,
                                   output_dir=str(folder / f"seed_{seed}"), sig_figs=12, show_console=False)
            names = list(fit.column_names)
            arr = np.asarray(fit.draws()).reshape(-1, len(names))
            lw = arr[:, names.index("lp__")] - arr[:, names.index("lp_approx__")]
            finite = np.isfinite(lw)
            # Rows are path-major (paths x draws).  Each path's mean log weight is its own
            # ELBO, a valid lower bound on log Z; the best path is the tightest bound.
            run["path_elbos"] = [float(np.mean(chunk[np.isfinite(chunk)])) if np.isfinite(chunk).any() else None
                                 for chunk in np.array_split(lw, paths)]
            if np.sum(finite) < 20:
                raise ValueError("Fewer than 20 finite importance weights")
            total += len(lw)
            lw = lw[finite]
            _, k = az.psislw(lw)
            w = np.exp(lw - logsumexp(lw))
            run.update(logz=float(logsumexp(lw) - np.log(len(lw))), elbo=float(lw.mean()),
                       ess=float(1 / np.sum(w**2)), pareto_k=float(k) if np.isfinite(k) else None,
                       draws=len(lw), discarded=int(np.sum(~finite)))
            for name, values in fit.stan_variables().items():
                if _kept(name):
                    retained[name].append(np.asarray(values)[finite])
            logs.append(lw)
        except Exception as ex:
            run["error"] = str(ex)
        runs.append(run)
    result = {"runs": runs, "seconds": time.monotonic() - t0, "method": "pathfinder_importance_sampling"}
    if not logs:
        return dict(result, status="failed", reliable=False), {}, None
    lw = np.concatenate(logs)
    np.savez_compressed(folder / "importance_weights.npz", **{f"start_{i}": v for i, v in enumerate(logs)})
    weights = np.exp(lw - logsumexp(lw))
    _, k = az.psislw(lw)
    ess = float(1 / np.sum(weights**2))
    reliable = (len(logs) == len(seeds) and ess >= minimum_ess and np.isfinite(k)
                and float(k) <= maximum_k and all(r.get("discarded", 1) == 0 for r in runs))
    result.update(status="ok", reliable=bool(reliable), ess=ess,
                  pareto_k=float(k) if np.isfinite(k) else None,
                  logz=float(logsumexp(lw) - np.log(total)), elbo=float(lw.mean()),
                  elbo_max=max(e for r in runs for e in r.get("path_elbos", []) if e is not None),
                  elbo_finite_draws_only=any(r.get("discarded", 0) for r in runs))
    sv = {name: np.concatenate(values) for name, values in retained.items()}
    np.savez_compressed(folder / "posterior_draws.npz", log_weight=lw, weight=weights,
                        **{k: v for k, v in sv.items() if k not in PREDICTIONS})
    if not keep_csv:
        for path in folder.glob("seed_*/*.csv"):
            path.unlink()
    return result, sv, weights


def nuts(model, data, inits, seed, folder, chains=4, warmup=1000, samples=1000,
         parallel_chains=4, adapt_delta=0.95):
    """NUTS reference fit.  Returns (result, draws, uniform weights).

    reliable = max R-hat <= 1.01, min bulk ESS >= 400, no divergences and no
    max-treedepth hits, over times, fractions and every Ne hyper/raw parameter.
    """
    result = dict(method="mcmc", sampling_settings=dict(warmup=warmup, samples=samples, chains=chains,
                                                        parallel_chains=parallel_chains, adapt_delta=adapt_delta))
    start = time.monotonic()
    try:
        fit = model.sample(data=data, inits=inits, seed=seed, chains=chains, parallel_chains=parallel_chains,
                           iter_warmup=warmup, iter_sampling=samples, adapt_delta=adapt_delta,
                           output_dir=str(folder / "chains"), show_progress=False)
        diagnostics = fit.method_variables()
        table = fit.summary()
        keys = [k for k in table.index
                if re.match(r"^(times\[|admixture_fractions\[|effective_N|mu_log|sigma_log|tau|Ne_raw)", k)]
        sub = table.loc[keys]
        rhat = float(sub["R_hat"].max()); ess = float(sub["ESS_bulk"].min())
        divergent = int(diagnostics["divergent__"].sum())
        depth_hits = int((diagnostics["treedepth__"] >= 10).sum())
        result.update(status="ok", reliable=bool(np.isfinite(rhat) and rhat <= 1.01 and ess >= 400
                                                 and divergent == 0 and depth_hits == 0),
                      max_rhat=rhat if np.isfinite(rhat) else None,
                      min_ess_bulk=ess if np.isfinite(ess) else None,
                      divergences=divergent, max_treedepth_hits=depth_hits)
        sv = {k: v for k, v in fit.stan_variables().items() if _kept(k)}
        weights = np.full(len(sv["cumulative_times"]), 1 / len(sv["cumulative_times"]))
    except Exception as ex:
        result.update(status="failed", reliable=False, error=str(ex))
        sv, weights = {}, None
    result["seconds"] = time.monotonic() - start
    return result, sv, weights
