#!/usr/bin/env python3
"""T_true against the two graphs that differ from it only in where b's minor source joins
the c-d clade (G1: with d; G2: with the ancestor of c and d), over many subsamples.

    python -B close_graphs.py run       # fits (resumable)
    python -B close_graphs.py report    # runs/default/close_graphs/report.md

Growth design, 100 cM (20 subsamples) and 1000 cM (5).  Per fit: best-path ELBO, the
importance-sampled log evidence (logZ) over all draws, and posterior means of the event
heights and f.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import json
import sys
import tempfile

import numpy as np
from scipy.special import logsumexp

import run
from methods.demography import DemographicTopology

OUT = run.HERE / "runs" / "default" / "close_graphs"
CASES = [("growth", 100, r) for r in range(20)] + [("growth", 1000, r) for r in range(5)]


def graphs():
    """Readable branch names; the events, in listing order, are named in EVENTS."""
    t = DemographicTopology(list("abcd"))
    t.add_admixture_event("b", "b1", "b2"); t.add_merge_event("a", "b1", "ab")
    t.add_merge_event("b2", "c", "cb2"); t.add_merge_event("cb2", "d", "cdb2"); t.add_merge_event("ab", "cdb2", "root")
    g1 = DemographicTopology(list("abcd"))
    g1.add_admixture_event("b", "b1", "b2"); g1.add_merge_event("a", "b1", "ab")
    g1.add_merge_event("b2", "d", "db2"); g1.add_merge_event("db2", "c", "cdb2"); g1.add_merge_event("ab", "cdb2", "root")
    g2 = DemographicTopology(list("abcd"))
    g2.add_admixture_event("b", "b1", "b2"); g2.add_merge_event("a", "b1", "ab")
    g2.add_merge_event("c", "d", "cd"); g2.add_merge_event("cd", "b2", "cdb2"); g2.add_merge_event("ab", "cdb2", "root")
    return {"T_true": t, "G1": g1, "G2": g2}


EVENTS = {"T_true": ["b admix", "a+b1", "b2+c", "+d", "root"],
          "G1": ["b admix", "a+b1", "b2+d", "+c", "root"],
          "G2": ["b admix", "a+b1", "c+d", "+b2", "root"]}


def fit(model, dem, data_args, seeds=run.SEEDS):
    import topology as tp
    from methods import fitting
    data = tp.stan_data(dem, *data_args)
    elbos, lws, heights, fs, weights = [], [], [], [], []
    for seed in seeds:
        with tempfile.TemporaryDirectory() as tmp:
            fit = model.pathfinder(data=data, inits=fitting.initial(data, tp.SPEC, seed), seed=seed, num_paths=1,
                                   draws=1000, num_single_draws=1000, psis_resample=False, calculate_lp=True,
                                   output_dir=tmp, sig_figs=12, show_console=False)
            cols = list(fit.column_names)
            arr = np.asarray(fit.draws()).reshape(-1, len(cols))
            lw = arr[:, cols.index("lp__")] - arr[:, cols.index("lp_approx__")]
            ok = np.isfinite(lw)
            sv = fit.stan_variables()
            elbos.append(float(lw[ok].mean())); lws.append(lw[ok])
            heights.append(np.asarray(sv["cumulative_times"])[ok]); fs.append(np.asarray(sv["admixture_fractions"])[ok])
    lw = np.concatenate(lws)
    w = np.exp(lw - logsumexp(lw))                          # importance weights over the pooled draws
    return dict(elbo=max(elbos), elbo_starts=elbos, logz=float(logsumexp(lw) - np.log(len(lw))),
                ess=float(1 / np.sum(w ** 2)),
                heights=(w[:, None] * np.concatenate(heights)).sum(0).tolist(),
                f=float((w * np.concatenate(fs).ravel()).sum()))


def run_fits(workers=10):
    from cmdstanpy import CmdStanModel
    import topology as tp
    model = CmdStanModel(stan_file=str(tp.STAN_FILE))
    todo = [(ds, g) for ds in CASES for g in graphs() if not (OUT / run.name(ds) / f"{g}.json").exists()]
    cache = {}

    def one(item):
        ds, g = item
        if ds not in cache:
            cache[ds] = run.observations(ds)
        r = fit(model, graphs()[g], cache[ds])
        dest = OUT / run.name(ds); dest.mkdir(parents=True, exist_ok=True)
        (dest / f"{g}.json").write_text(json.dumps(r))
    for ds in CASES:                                         # load data once, serially
        cache[ds] = run.observations(ds)
    with ThreadPoolExecutor(workers) as pool:
        list(pool.map(one, todo))
    print(f"{len(todo)} fits done")


def report():
    rows = {(ds, g): json.loads((OUT / run.name(ds) / f"{g}.json").read_text()) for ds in CASES for g in graphs()}
    lines = ["# T_true against G1 and G2 (where b's minor source joins the c-d clade)", "",
             "Growth design.  ELBO = best path of 3 starts; logZ = importance-sampled log evidence over all 3000 draws "
             "(an estimate of log p(x | G), at most the true value in expectation).  Heights and f: importance-weighted "
             "posterior means.  Truth: b admix 20 (f = 0.7 to the a side), a+b1 60, b2+c 100, +d 200, root 400.", ""]
    for cm in sorted({ds[1] for ds in CASES}):
        cases = [ds for ds in CASES if ds[1] == cm]
        lines += [f"## {cm} cM ({len(cases)} subsamples)", "",
                  "| | G2 − T_true | G1 − T_true |", "|---|---|---|"]
        for stat in ("elbo", "logz"):
            cells = []
            for g in ("G2", "G1"):
                d = np.array([rows[(ds, g)][stat] - rows[(ds, "T_true")][stat] for ds in cases])
                cells.append(f"median {np.median(d):+.1f} (range {d.min():+.1f} to {d.max():+.1f}); "
                             f"above T_true in {np.mean(d > 0):.0%}")
            lines.append(f"| {stat} | " + " | ".join(cells) + " |")
        lines += ["", "| subsample | graph | ELBO | logZ | ESS | event heights (posterior mean) | f |",
                  "|---:|---|---:|---:|---:|---|---:|"]
        for ds in cases:
            for g in graphs():
                r = rows[(ds, g)]
                h = ", ".join(f"{e} {x:.0f}" for e, x in zip(EVENTS[g], r["heights"]))
                lines.append(f"| {ds[2]} | {g} | {r['elbo']:.1f} | {r['logz']:.1f} | {r['ess']:.0f} | {h} | {r['f']:.2f} |")
        lines.append("")
    (OUT / "report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    {"run": run_fits, "report": report}[sys.argv[1]]()
