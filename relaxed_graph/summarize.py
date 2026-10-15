"""Report and figures for runs/<out>/fits.

  report.md             per genome length: exact EM (theta_hat, graph posterior), profile gaps,
                        relaxed MC-EM per temperature (theta_hat, read-off, distance to exact,
                        exact objective at the relaxed theta_hat, convergence)
  graph_posteriors.png  p(G | x, theta_hat): exact vs relaxed read-off per temperature
  em_traces.png         relaxed MC-EM: mean log-likelihood and theta by iteration
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import data as dt
from models import compile_model, at

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
EXACT = "#3d3d3a"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "axes.axisbelow": True})


def tag(key):
    lam, start = key
    return f"relaxed λ={lam:g}" + ("" if start == "neutral" else f", from exact θ̂")


def lam_color(key, keys):
    """Sequential blue by temperature: darker = colder (smaller lambda)."""
    cmap = plt.get_cmap("Blues")
    order = sorted({k[0] for k in keys}, reverse=True)
    return cmap(0.35 + 0.6 * order.index(key[0]) / max(1, len(order) - 1))


def graph_probs(sk, vertex_values):
    return {k: float(np.sum(np.asarray(vertex_values)[ix])) for k, ix in sk.graphs().items()}


def load(dest):
    out = {}
    for name in ("discrete", "profile"):
        if (dest / f"{name}.json").exists():
            out[name] = json.loads((dest / f"{name}.json").read_text())
    out["relaxed"] = {}
    for f in sorted(dest.glob("relaxed_lambda_*.json")):
        r = json.loads(f.read_text())
        out["relaxed"][(r["temperature"], r.get("start", "neutral"))] = r
    return out


def times(theta):
    return np.round(np.cumsum(theta["times"])).astype(int).tolist()


def build(c, sk, out):
    dest = out / "summary"; dest.mkdir(parents=True, exist_ok=True)
    discrete_model = compile_model("discrete")
    lines = ["# Relaxed graph vs exact EM (3 leaves, skeleton [admixture, merge, merge, root])", "",
             f"Data: `{c['trees_from']}` (true IBD, {len(c['bin_edges_cm']) - 1} bins from "
             f"{c['bin_edges_cm'][0]} cM, SNP {'on' if c['use_snp'] else 'off'}); theta = (event times, f, N) "
             f"as constants.  {len(sk.vertices)} vertices = {len(sk.graphs())} graphs.  Truth (outside the "
             "skeleton): loop on b at 10-25, admixture of b at 30 (0.75 from the a side), a+b1 at 50, "
             "b2+c at 100, root 300, N 15,000.", ""]
    panels = []
    for cm in c["genome_cm"]:
        folder = out / "fits" / f"cm_{int(cm):04d}" / "rep_000"
        if not folder.exists():
            continue
        r = load(folder)
        if "discrete" not in r:
            continue
        exact = graph_probs(sk, r["discrete"]["vertex_posterior"])
        th = r["discrete"]["theta"]
        lines += [f"## {cm} cM", "",
                  f"**Exact EM**: log p(x | theta_hat) = {r['discrete']['log_marginal']:.2f}; event times "
                  f"{times(th)}, f {th['f']:.3f}, N {th['N']:,.0f} ({r['discrete']['seconds']:.0f} s).", ""]
        obs, _ = dt.subsample(c, out, cm, 0)
        base = dt.stan_data(c, sk, obs)
        vdata = dict(base, V=len(sk.vertices), vertex_admix=[v[0] + 1 for v in sk.vertices],
                     vertex_merge=[[p + 1 for p in v[1:]] for v in sk.vertices],
                     log_prior=[-np.log(len(sk.vertices))] * len(sk.vertices))
        profile = {}
        if "profile" in r:
            ll = np.asarray(r["profile"]["vertex_loglik"])
            profile = {k: float(ll[ix].max()) for k, ix in sk.graphs().items()}
        best_profile = max(profile.values()) if profile else None
        relaxed = {}
        for lam, rr in sorted(r["relaxed"].items(), key=lambda kv: (kv[0][1] != "neutral", -kv[0][0])):
            probs = graph_probs(sk, rr["vertex_posterior"])
            tv = 0.5 * sum(abs(probs[k] - exact[k]) for k in exact)
            # the exact objective log p(x | theta) at the relaxed estimate
            ex = at(discrete_model, vdata, rr["theta"])
            lpx = float(np.log(np.sum(np.exp(ex["loglik"] - ex["loglik"].max()))) + ex["loglik"].max()
                        - np.log(len(sk.vertices)))
            relaxed[lam] = dict(probs=probs, tv=tv, lpx=lpx, theta=rr["theta"], iterations=len(rr["trace"]),
                                trace=rr["trace"], rhat=rr["max_rhat"], seconds=rr["seconds"])
        names = sorted(exact, key=lambda k: -exact[k])[:4]
        for lam in relaxed:
            names += [k for k in sorted(relaxed[lam]["probs"], key=lambda k: -relaxed[lam]["probs"][k])[:2]
                      if k not in names]
        head = "| graph | exact EM | " + " | ".join(tag(lam) for lam in relaxed) + \
               (" | profile Δℓ |" if profile else " |")
        lines += [head, "|---|---:|" + "---:|" * len(relaxed) + ("---:|" if profile else "")]
        for k in names:
            cells = [f"{exact[k]:.3f}"] + [f"{relaxed[lam]['probs'][k]:.3f}" for lam in relaxed]
            if profile:
                cells.append(f"{profile[k] - best_profile:+.1f}")
            lines.append(f"| {k} | " + " | ".join(cells) + " |")
        if relaxed:
            lines += ["", "| λ | event times | f | N | exact log p(x \\| θ) at this θ | TV to exact | "
                          "EM iterations | max R-hat | seconds |", "|---:|---|---:|---:|---:|---:|---:|---:|---:|"]
            for lam, v in relaxed.items():
                lines.append(f"| {tag(lam)} | {times(v['theta'])} | {v['theta']['f']:.3f} | {v['theta']['N']:,.0f} | "
                             f"{v['lpx']:.2f} | {v['tv']:.3f} | {v['iterations']} | {v['rhat']:.3f} | "
                             f"{v['seconds']:.0f} |")
        lines.append("")
        panels.append((cm, exact, relaxed, names))
    (dest / "report.md").write_text("\n".join(lines) + "\n")
    if panels:
        figures(panels, dest)
    print(f"Summary in {dest}")


def figures(panels, dest):
    lams = sorted({lam for _, _, rel, _ in panels for lam in rel})
    fig, axes = plt.subplots(len(panels), 1, figsize=(13, 3.6 * len(panels)), squeeze=False)
    for ax, (cm, exact, relaxed, names) in zip(axes[:, 0], panels):
        x = np.arange(len(names)); n = 1 + len(relaxed); width = 0.8 / n
        ax.bar(x - 0.4 + width / 2, [exact[k] for k in names], width * 0.9, color=EXACT, label="exact EM")
        for i, lam in enumerate(relaxed):
            ax.bar(x - 0.4 + width * (i + 1.5), [relaxed[lam]["probs"][k] for k in names], width * 0.9,
                   color=lam_color(lam, lams), label=tag(lam))
        ax.set_xticks(x, [k.replace(" ; ", "\n") for k in names], fontsize=7)
        ax.set_ylim(0, 1.05); ax.set_ylabel("p(G | x, θ̂)")
        ax.set_title(f"{cm} cM", loc="left", fontsize=10)
        ax.legend(frameon=False, fontsize=8, loc="upper right", ncol=min(5, n))
    fig.suptitle("Graph posterior at the EM estimate: exact (sum over the 108 vertices) vs relaxed read-off",
                 x=0.01, ha="left", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(dest / "graph_posteriors.png", dpi=180, bbox_inches="tight"); plt.close(fig)
    fig, axes = plt.subplots(len(panels), 3, figsize=(15, 3.2 * len(panels)), squeeze=False)
    for row, (cm, _, relaxed, _) in zip(axes, panels):
        for lam, v in relaxed.items():
            it = [t["iteration"] for t in v["trace"]]
            kw = dict(color=lam_color(lam, lams), lw=2, marker="o", ms=3, label=tag(lam),
                      ls="-" if lam[1] == "neutral" else "--")
            row[0].plot(it, [t["mean_loglik"] for t in v["trace"]], **kw)
            row[1].plot(it, [t["theta"]["f"] for t in v["trace"]], **kw)
            row[2].plot(it, [np.cumsum(t["theta"]["times"])[0] for t in v["trace"]], **kw)
        row[0].set_ylabel("mean ℓ over E-step draws"); row[1].set_ylabel("f"); row[2].set_ylabel("admixture time")
        row[0].set_title(f"{cm} cM", loc="left", fontsize=10)
        for ax in row:
            ax.set_xlabel("EM iteration")
        if relaxed:
            row[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Relaxed Monte Carlo EM traces (θ entering each E-step)", x=0.01, ha="left", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(dest / "em_traces.png", dpi=180, bbox_inches="tight"); plt.close(fig)
