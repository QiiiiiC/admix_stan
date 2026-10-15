#!/usr/bin/env python3
"""Stages: extract, discrete, profile, relaxed, summarize.  Every result is one JSON file,
skipped once written, so stages resume and can run in parallel.

  discrete   exact EM: theta_hat = argmax log sum_v p(v) exp l(theta, v); p(v | x, theta_hat)
  profile    every vertex's own maximum likelihood max_theta l(theta, v)
  relaxed    Monte Carlo EM on the Gumbel-softmax relaxation at each temperature lambda
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import time

import numpy as np

import data as dt
from models import compile_model
from skeleton import Skeleton
from methods.utils import save_json

logging.getLogger("cmdstanpy").setLevel(logging.ERROR)
HERE = Path(__file__).resolve().parent


def folder(out, cm, rep):
    return out / "fits" / f"cm_{int(cm):04d}" / f"rep_{rep:03d}"


def random_theta(rng, K):
    """Starting values: increments log-uniform on [5, 300] generations, f uniform, N
    log-uniform on [3,000, 60,000]."""
    return dict(times=np.exp(rng.uniform(np.log(5), np.log(300), K)).tolist(), f=float(rng.uniform(0.1, 0.9)),
                N=float(np.exp(rng.uniform(np.log(3000), np.log(60000)))))


def theta_of(fit):
    v = fit.stan_variables()
    return dict(times=np.atleast_1d(v["times"]).tolist(), f=float(v["f"]), N=float(v["N"]))


def optimize(model, data, inits, seed):
    """Best of several L-BFGS runs (maximum likelihood: no Jacobian)."""
    best = None
    for k, init in enumerate(inits):
        try:
            fit = model.optimize(data=data, inits=init, seed=seed + k, algorithm="lbfgs", jacobian=False,
                                 iter=5000, show_console=False, sig_figs=17)
        except RuntimeError:
            continue
        lp = float(fit.optimized_params_dict["lp__"])
        if best is None or lp > best[0]:
            best = (lp, fit)
    return best


def vertex_data(base, sk, vertices):
    return dict(base, V=len(vertices), vertex_admix=[v[0] + 1 for v in vertices],
                vertex_merge=[[p + 1 for p in v[1:]] for v in vertices],
                log_prior=[-np.log(len(sk.vertices))] * len(vertices))


def graph_table(sk, vertex_values, combine):
    """Per distinct graph: combine the values of its vertices (sum for probabilities,
    max for profile log-likelihoods)."""
    return {sk.label(sk.vertices[ix[0]]): combine([vertex_values[i] for i in ix]) for ix in sk.graphs().values()}


# ---------------------------------------------------------------------------
def discrete(c, sk, base, dest, rng):
    t0 = time.monotonic()
    model = compile_model("discrete")
    inits = [random_theta(rng, sk.K) for _ in range(c["optimize_starts"])]
    lp, fit = optimize(model, vertex_data(base, sk, sk.vertices), inits, c["seed"])
    post = np.atleast_1d(fit.stan_variables()["posterior"])
    ll = np.atleast_1d(fit.stan_variables()["loglik"])
    save_json(dest / "discrete.json", dict(
        method="discrete EM", log_marginal=lp, theta=theta_of(fit), vertex_posterior=post.tolist(),
        vertex_loglik=ll.tolist(), graph_posterior=graph_table(sk, post, sum), seconds=time.monotonic() - t0))


def profile(c, sk, base, dest, rng):
    t0 = time.monotonic()
    model = compile_model("discrete")
    ll, thetas = [], []
    for v in sk.vertices:
        inits = [random_theta(rng, sk.K) for _ in range(3)]
        lp, fit = optimize(model, vertex_data(base, sk, [v]), inits, c["seed"])
        ll.append(lp + np.log(len(sk.vertices)))             # undo the log prior: this is l(theta_v, v)
        thetas.append(theta_of(fit))
    save_json(dest / "profile.json", dict(
        method="profile", vertex_loglik=ll, vertex_theta=thetas, graph_loglik=graph_table(sk, ll, max),
        seconds=time.monotonic() - t0))


def result_name(lam, start):
    return f"relaxed_lambda_{lam:g}" + ("" if start == "neutral" else f"_start_{start}") + ".json"


def relaxed(c, sk, base, dest, lam, rng, start="neutral"):
    """start: "neutral" (the same generic theta for every lambda) or "exact" (the exact
    EM's theta_hat: a diagnostic of whether the relaxed EM keeps the right answer)."""
    t0 = time.monotonic()
    e, m = c["em"], sk.n_merge
    estep, mstep = compile_model("relaxed_estep"), compile_model("relaxed_mstep")
    if start == "exact":
        theta = json.loads((dest / "discrete.json").read_text())["theta"]
    else:
        theta = dict(times=[10.0, 40.0, 100.0, 200.0][:sk.K], f=0.5, N=10000.0)
    trace = []
    for it in range(e["iterations"]):
        fit = estep.sample(data=dict(base, **{"lambda": lam}, **theta), chains=e["chains"], parallel_chains=e["chains"],
                           iter_warmup=e["warmup"], iter_sampling=e["samples"], seed=c["seed"] + it,
                           show_progress=False, show_console=False, adapt_delta=0.9)
        sv = fit.stan_variables()
        n = len(sv["loglik"])
        pick = rng.choice(n, size=min(e["mstep_draws"], n), replace=False)
        wa = sv["w_admix"][pick]; wm = sv["w_merge"][pick]
        summary = fit.summary()
        trace.append(dict(iteration=it, theta=theta, mean_loglik=float(np.mean(sv["loglik"])),
                          max_rhat=float(np.nanmax(summary["R_hat"].values)),
                          divergences=int(fit.divergences.sum()) if fit.divergences is not None else None))
        lp, mfit = optimize(mstep, dict(base, S=len(pick), w_admix=wa.tolist(), w_merge=wm.tolist()), [theta],
                            c["seed"])
        new = theta_of(mfit)
        change = max(abs(np.log(np.r_[new["times"], new["f"], new["N"]]) -
                         np.log(np.r_[theta["times"], theta["f"], theta["N"]])))
        theta = new
        print(f"lambda={lam} iteration {it}: mean loglik {trace[-1]['mean_loglik']:.1f}, "
              f"max |dlog theta| {change:.4f}", flush=True)
        if change < e["tolerance"]:
            break
    # final E-step at theta_hat: read-off of the hard graph
    fit = estep.sample(data=dict(base, **{"lambda": lam}, **theta), chains=e["chains"], parallel_chains=e["chains"],
                       iter_warmup=e["warmup"], iter_sampling=e["samples"], seed=c["seed"] + 999,
                       show_progress=False, show_console=False, adapt_delta=0.9)
    sv = fit.stan_variables()
    hard = [(int(a) - 1, *[int(p) - 1 for p in ps]) for a, ps in zip(sv["hard_admix"], sv["hard_merge"])]
    index = {v: i for i, v in enumerate(sk.vertices)}
    counts = np.bincount([index[h] for h in hard], minlength=len(sk.vertices)) / len(hard)
    chain = np.repeat(np.arange(e["chains"]), e["samples"])
    per_chain = {int(k): graph_table(sk, np.bincount([index[h] for h, ch in zip(hard, chain) if ch == k],
                                                     minlength=len(sk.vertices)) / e["samples"], sum)
                 for k in range(e["chains"])}
    save_json(dest / result_name(lam, start), dict(
        method="relaxed MC-EM", temperature=lam, start=start, theta=theta, trace=trace, vertex_posterior=counts.tolist(),
        graph_posterior=graph_table(sk, counts, sum), graph_posterior_by_chain=per_chain,
        soft_weight_mean=dict(admix=sv["w_admix"].mean(0).tolist(), merge=sv["w_merge"].mean(0).tolist()),
        max_rhat=float(np.nanmax(fit.summary()["R_hat"].values)), seconds=time.monotonic() - t0))


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["extract", "discrete", "profile", "relaxed", "summarize"])
    p.add_argument("--out", type=Path, default=HERE / "runs" / "default")
    p.add_argument("--cm", type=float, nargs="+")
    p.add_argument("--lambdas", type=float, nargs="+")
    p.add_argument("--start", choices=["neutral", "exact"], default="neutral")
    a = p.parse_args(); out = a.out.resolve()
    c = dt.load_config(); sk = Skeleton(c["leaves"])
    if a.stage == "extract":
        dt.extract(c, out); return
    if a.stage == "summarize":
        from summarize import build
        build(c, sk, out); return
    for cm in a.cm or c["genome_cm"]:
        for rep in range(c["replicates"]):
            dest = folder(out, cm, rep); dest.mkdir(parents=True, exist_ok=True)
            obs, chosen = dt.subsample(c, out, cm, rep)
            base = dt.stan_data(c, sk, obs)
            save_json(dest / "subsample.json", dict(genome_cm=cm, replicate=rep, blocks=chosen))
            rng = np.random.default_rng([c["seed"], int(cm), rep, ["discrete", "profile", "relaxed"].index(a.stage)])
            if a.stage == "discrete" and not (dest / "discrete.json").exists():
                discrete(c, sk, base, dest, rng)
            elif a.stage == "profile" and not (dest / "profile.json").exists():
                profile(c, sk, base, dest, rng)
            elif a.stage == "relaxed":
                for lam in a.lambdas or c["lambdas"]:
                    if not (dest / result_name(lam, a.start)).exists():
                        relaxed(c, sk, base, dest, lam, rng, a.start)


if __name__ == "__main__":
    main()
