#!/usr/bin/env python3
"""Search the graph space on the four-leaf data, against the full enumeration.

    python -B run.py fit       # Pathfinder ELBO of every graph in the space, per dataset (resumable)
    python -B run.py search    # walks and greedy searches replayed on the cached ELBOs
    python -B run.py report    # runs/default/report.md and figures

The space is every graph on a, b, c, d with at most one admixture, no bubble and no
ghost source (153 graphs).  Each graph is fitted once with the per-branch model
(branch_model/), Pathfinder 3 starts x 1 path; its score is the best-path ELBO.

A search only ever needs the score of graphs it proposes, and a graph's score is fixed
once fitted, so a search replayed on the enumeration's cache is exactly the search that
would have fitted on demand; "fits" counts the distinct graphs it scored.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import logging
from pathlib import Path
import sys
import tempfile
import time

import numpy as np
from scipy.special import logsumexp

import graphs as gr

HERE = Path(__file__).resolve().parent
STUDY = gr.ROOT / "simulation_study" / "four_leaf_b_admixture_growth_ne"
DESIGNS = {"growth": STUDY, "varying": STUDY.parent / "four_leaf_b_admixture_varying_ne"}
DATASETS = [(d, cm, 0) for d in DESIGNS for cm in (100, 1000)]
SPACE = gr.Space(max_admix=1)
LEAVES = "abcd"
SEEDS = (1, 7, 13)
logging.getLogger("cmdstanpy").setLevel(logging.ERROR)


def name(ds):
    return f"{ds[0]}_cm{ds[1]:04d}_rep{ds[2]:03d}"


def file_of(key):
    return hashlib.sha1(key.encode()).hexdigest()[:16] + ".json"


# ---------------------------------------------------------------------------
def observations(ds):
    sys.path.insert(0, str(STUDY))
    import study as s
    from methods.blocks import load_blocks
    design, cm, rep = ds
    c = s.load_config(DESIGNS[design] / "config.json")
    blocks = load_blocks(DESIGNS[design] / "runs" / "default" / "pool_000", c["blocks"], ["true"])
    obs = s.aggregate(blocks, s.selections(c, 0, cm)[rep], c, "true")
    return obs, s.edges(c), s.n_haploid(c), c["T_max"]


def fit_graph(model, g, obs, edges, n_haploid, T_max):
    import topology as tp
    from methods import fitting
    data = tp.stan_data(gr.to_topology(g), obs, edges, n_haploid, T_max)
    starts, t0 = [], time.monotonic()
    for seed in SEEDS:
        with tempfile.TemporaryDirectory() as tmp:
            try:
                fit = model.pathfinder(data=data, inits=fitting.initial(data, tp.SPEC, seed), seed=seed, num_paths=1,
                                       draws=1000, num_single_draws=1000, psis_resample=False, calculate_lp=True,
                                       output_dir=tmp, sig_figs=12, show_console=False)
                cols = list(fit.column_names)
                arr = np.asarray(fit.draws()).reshape(-1, len(cols))
                lw = arr[:, cols.index("lp__")] - arr[:, cols.index("lp_approx__")]
                lw = lw[np.isfinite(lw)]
                starts.append(dict(seed=seed, elbo=float(lw.mean()), logz=float(logsumexp(lw) - np.log(len(lw)))))
            except Exception as ex:
                starts.append(dict(seed=seed, error=str(ex)[-300:]))
    ok = [s for s in starts if "elbo" in s]
    return dict(key=g.key(), n_admix=g.n_admix(), starts=starts, seconds=time.monotonic() - t0,
                elbo=max(s["elbo"] for s in ok) if ok else None)


def fit(out, workers):
    from cmdstanpy import CmdStanModel
    import topology as tp
    model = CmdStanModel(stan_file=str(tp.STAN_FILE))
    seen, _ = gr.enumerate_space(SPACE, LEAVES)
    for ds in DATASETS:
        dest = out / name(ds) / "fits"; dest.mkdir(parents=True, exist_ok=True)
        obs, edges, n_haploid, T_max = observations(ds)
        todo = [g for k, g in seen.items() if not (dest / file_of(k)).exists()]
        t0 = time.monotonic()

        def one(g):
            r = fit_graph(model, g, obs, edges, n_haploid, T_max)
            (dest / file_of(r["key"])).write_text(json.dumps(r))
        with ThreadPoolExecutor(workers) as pool:
            list(pool.map(one, todo))
        print(f"{name(ds)}: {len(todo)} graphs fitted in {time.monotonic() - t0:.0f} s", flush=True)


# ---------------------------------------------------------------------------
def load_scores(out, ds):
    rows = [json.loads(f.read_text()) for f in (out / name(ds) / "fits").glob("*.json")]
    return {r["key"]: (r["elbo"] if r["elbo"] is not None else -np.inf) for r in rows}, rows


def log_prior(seen):
    """Uniform over the number of admixtures, then uniform within it."""
    count = {}
    for g in seen.values():
        count[g.n_admix()] = count.get(g.n_admix(), 0) + 1
    return {k: -np.log(len(count)) - np.log(count[g.n_admix()]) for k, g in seen.items()}


def walk(nbrs, score, prior, start, beta, steps, rng):
    """Metropolis-Hastings over graphs on cached scores.  Returns the order in which
    graphs were first scored (the fits a live run would make)."""
    total = {k: sum(m for m, _ in n.values()) for k, n in nbrs.items()}
    cur, scored = start, [start]
    known = {start}
    for _ in range(steps):
        keys = list(nbrs[cur]); mult = np.array([nbrs[cur][k][0] for k in keys], float)
        new = keys[rng.choice(len(keys), p=mult / mult.sum())]
        if new not in known:
            known.add(new); scored.append(new)
        log_a = (beta * (score[new] - score[cur]) + prior[new] - prior[cur]
                 + np.log(nbrs[new][cur][0] / total[new]) - np.log(nbrs[cur][new][0] / total[cur]))
        if np.log(rng.uniform()) < log_a:
            cur = new
    return scored


def greedy(nbrs, score, start, allowed=lambda k: True):
    """Score every neighbour, move to the best if it improves; stop at a local optimum."""
    cur, scored, known = start, [start], {start}
    while True:
        cand = [k for k in nbrs[cur] if allowed(k)]
        for k in cand:
            if k not in known:
                known.add(k); scored.append(k)
        best = max(cand, key=score.get)
        if score[best] <= score[cur]:
            return scored
        cur = best


def search(out, chains, steps):
    seen, nbrs = gr.enumerate_space(SPACE, LEAVES)
    _, nbrs_rel = gr.enumerate_space(SPACE, LEAVES, relocate=True)
    trees = [k for k, g in seen.items() if g.n_admix() == 0]
    prior = log_prior(seen)
    for ds in DATASETS:
        score, _ = load_scores(out, ds)
        assert set(score) == set(seen), "fit stage incomplete"
        rng = np.random.default_rng([1, ds[1]])
        runs = {}
        for beta in (1.0, 0.3, 0.1):
            runs[f"walk, beta={beta:g}"] = [walk(nbrs, score, prior, trees[rng.integers(len(trees))], beta, steps, rng)
                                           for _ in range(chains)]
        for beta in (1.0, 0.3):
            runs[f"walk + relocate, beta={beta:g}"] = [
                walk(nbrs_rel, score, prior, trees[rng.integers(len(trees))], beta, steps, rng) for _ in range(chains)]
        runs["greedy"] = [greedy(nbrs, score, t) for t in trees]
        runs["greedy + relocate"] = [greedy(nbrs_rel, score, t) for t in trees]
        stepwise = []
        for t in trees:                      # TreeMix-like: best tree first, then the full space
            first = greedy(nbrs, score, t, allowed=lambda k: seen[k].n_admix() == 0)
            best_tree = max(first, key=score.get)
            rest = greedy(nbrs, score, best_tree)
            stepwise.append(first + [k for k in rest if k not in first])
        runs["stepwise (trees, then all)"] = stepwise
        (out / name(ds) / "search.json").write_text(json.dumps(runs))
        print(f"{name(ds)}: searches done", flush=True)


# ---------------------------------------------------------------------------
def report(out):
    from summarize import build
    seen, nbrs = gr.enumerate_space(SPACE, LEAVES)
    build(out, DATASETS, name, seen, log_prior(seen), load_scores, nbrs)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("stage", choices=["fit", "search", "report"])
    p.add_argument("--out", type=Path, default=HERE / "runs" / "default")
    p.add_argument("--workers", type=int, default=10)
    p.add_argument("--chains", type=int, default=200)
    p.add_argument("--steps", type=int, default=300)
    a = p.parse_args()
    out = a.out.resolve()
    if a.stage == "fit":
        fit(out, a.workers)
    elif a.stage == "search":
        search(out, a.chains, a.steps)
    else:
        report(out)


if __name__ == "__main__":
    main()
