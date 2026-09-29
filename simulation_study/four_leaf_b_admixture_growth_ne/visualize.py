"""Tables, report and figures.

  topology_rank        ΔELBO(T_true − T_null) (best-path ELBO) per model vs length, and the
                       rate at which T_true ranks first (Wilson 95%); one row per IBD source
  parameters_<src>     relative error of the 16 T_true parameters: Pathfinder for both models
                       (all subsamples) and, for the NUTS source, NUTS for the mixed model
  spectrum_<L>cM_<src> observed IBD counts in every bin and pair across all subsamples at the
                       longest length, with each model's fitted expected counts
  residuals_<L>cM_<src> Pearson residuals of those counts per bin and pair
  hapibd_detection     hap-IBD segment counts over true counts per bin, all blocks pooled
  elbo_vs_logz         exact log Z (NUTS + bridge sampling) against the best-path ELBO for the
                       mixed model on both graphs: the gap, ΔlogZ vs ΔELBO, identification
  topology_true        the generating graph with every branch's Ne
"""
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import norm

from study import EVENTS, GRAPHS, MODEL_NAMES, NODES, POPS, branch_ne, build, edges, n_haploid, parameters, save_json
from methods.ibd import pair_counts

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
# Reference categorical slots 1, 2, 3 (validated all-pairs)
COLORS = {"ibd": "#2a78d6", "mixed": "#eb6834", "nuts": "#1baf7a"}
LABELS = {"ibd": "IBD-only (Pathfinder)", "mixed": "Mixed (Pathfinder)", "nuts": "Mixed (NUTS)"}
SOURCE_LABELS = {"true": "true IBD", "hapibd": "hap-IBD"}
GRAPH_COLORS = {"T_true": "#eb6834", "T_null": "#2a78d6"}
BACK = dict(boxstyle="round,pad=0.15", fc="white", ec="none")
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "axes.axisbelow": True})


def csv_write(path, rows, fields):
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fields, extrasaction="ignore"); w.writeheader(); w.writerows(rows)


def save(fig, dest, name):
    Path(dest).mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(Path(dest) / f"{name}.{ext}", dpi=200 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def wilson(k, n, level=0.95):
    if n == 0:
        return np.nan, np.nan
    z = norm.ppf(0.5 + level / 2); p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half


def load(out, stage):
    res = {}
    for f in sorted((out / "pool_000" / stage).glob("cm_*/rep_*/*/*/*/result.json")):
        r = json.loads(f.read_text())
        res[(r["genome_cm"], r["replicate"], r["source"], r["model"], r["graph"])] = r
    return res


def build_summary(c, out):
    out = Path(out); dest = out / "summary"; fig_dir = dest / "figures"; fig_dir.mkdir(parents=True, exist_ok=True)
    pf, nuts = load(out, "pathfinder"), load(out, "nuts")
    fit_rows = [dict(stage="pathfinder", genome_cm=k[0], replicate=k[1], source=k[2], model=k[3], graph=k[4],
                     status=r["status"],
                     elbo=r.get("elbo"), elbo_max=r.get("elbo_max"), logz_is=r.get("logz"), ess=r.get("ess"),
                     pareto_k=r.get("pareto_k"), seconds=r["seconds"]) for k, r in pf.items()]
    fit_rows += [dict(stage="nuts", genome_cm=k[0], replicate=k[1], source=k[2], model=k[3], graph=k[4],
                      status=r["status"],
                      logz_bridge=r.get("evidence", {}).get("logz"), bridge_error=r.get("evidence", {}).get("error"),
                      max_rhat=r.get("max_rhat"), min_ess_bulk=r.get("min_ess_bulk"), divergences=r.get("divergences"),
                      max_treedepth_hits=r.get("max_treedepth_hits"), seconds=r["seconds"]) for k, r in nuts.items()]
    csv_write(dest / "fits.csv", fit_rows, ["stage", "genome_cm", "replicate", "source", "model", "graph", "status", "elbo",
                                            "elbo_max", "logz_is", "ess", "pareto_k", "logz_bridge", "bridge_error",
                                            "max_rhat", "min_ess_bulk", "divergences", "max_treedepth_hits", "seconds"])
    topo = topology_table(c, pf)
    csv_write(dest / "topology.csv", topo["rows"], ["genome_cm", "replicate", "source", "model", "true_elbo_max", "null_elbo_max",
                                                     "delta_elbo_max", "true_first"])
    params = parameter_rows(pf, nuts)
    csv_write(dest / "parameter_estimates.csv", params, ["source", "method", "genome_cm", "replicate", "parameter", "mean",
                                                         "lower", "upper", "truth", "relative_error"])
    ev = evidence_rows(pf, nuts)
    csv_write(dest / "elbo_vs_logz.csv", ev, ["genome_cm", "replicate", "source", "graph", "logz", "logz_error", "elbo_max",
                                              "elbo", "gap"])
    save_json(dest / "summary_status.json", dict(
        pathfinder=dict(fits=len(pf), failed=sum(r["status"] != "ok" for r in pf.values())),
        nuts=dict(fits=len(nuts), failed=sum(r["status"] != "ok" for r in nuts.values()),
                  converged=sum(bool(r.get("reliable")) for r in nuts.values()))))
    topology_figure_rank(c, topo, fig_dir)
    for source in c["sources"]:
        parameter_figure(c, [r for r in params if r["source"] == source], fig_dir, source)
        spectrum_figures(c, {k: r for k, r in pf.items() if k[2] == source}, fig_dir, source)
    detection = hapibd_detection(c, out) if "hapibd" in c["sources"] else None
    if detection:
        csv_write(dest / "hapibd_detection.csv", detection["rows"], ["bin_lower_cm", "bin_upper_cm", "pair", "true",
                                                                      "hapibd", "ratio"])
        detection_figure(c, detection, fig_dir)
    if ev:
        elbo_logz_figure(c, ev, fig_dir)
    true_graph_figure(c, fig_dir)
    report(c, dest, topo, params, ev, pf, nuts, detection)
    print(f"Summarised {len(pf)} Pathfinder and {len(nuts)} NUTS fits into {dest}")


# ---------------------------------------------------------------------------
def topology_table(c, pf):
    rows, table = [], []
    for (cm, rep, source, model, graph), r in pf.items():
        if graph != "T_true":
            continue
        n = pf.get((cm, rep, source, model, "T_null"))
        if r["status"] != "ok" or not n or n["status"] != "ok":
            continue
        d = r["elbo_max"] - n["elbo_max"]
        rows.append(dict(genome_cm=cm, replicate=rep, source=source, model=model, true_elbo_max=r["elbo_max"],
                         null_elbo_max=n["elbo_max"], delta_elbo_max=d, true_first=d > 0))
    for source in c["sources"]:
        for model in MODEL_NAMES:
            for cm in c["genome_cm"]:
                g = [x for x in rows if x["source"] == source and x["model"] == model and x["genome_cm"] == cm]
                if g:
                    k = sum(x["true_first"] for x in g); d = np.array([x["delta_elbo_max"] for x in g])
                    table.append(dict(source=source, model=model, genome_cm=cm, n=len(g), rate=k / len(g),
                                      ci=wilson(k, len(g)), median=float(np.median(d)),
                                      q25=float(np.percentile(d, 25)), q75=float(np.percentile(d, 75)),
                                      minimum=float(d.min())))
    return dict(rows=rows, table=table)


def topology_figure_rank(c, topo, dest):
    lengths = c["genome_cm"]; x = np.arange(len(lengths)); width = 0.36; sources = c["sources"]
    fig, axes = plt.subplots(len(sources), 2, figsize=(14, 4.4 * len(sources)), sharey="col", squeeze=False)
    for (a1, a2), source in zip(axes, sources):
        rows = [r for r in topo["rows"] if r["source"] == source]
        table = [r for r in topo["table"] if r["source"] == source]
        for m, model in enumerate(MODEL_NAMES):
            data = [[r["delta_elbo_max"] for r in rows if r["model"] == model and r["genome_cm"] == cm] for cm in lengths]
            pos = [k for k, d in enumerate(data) if d]
            if not pos:
                continue
            bp = a1.boxplot([data[k] for k in pos], positions=[k + (m - 0.5) * width for k in pos], widths=width * 0.85,
                            patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1.3),
                            whiskerprops=dict(color=MUTED), capprops=dict(color=MUTED))
            for box in bp["boxes"]:
                box.set(facecolor=COLORS[model], alpha=0.6, edgecolor=MUTED)
            t = [r for r in table if r["model"] == model]
            xs = [lengths.index(r["genome_cm"]) for r in t]
            y = np.array([100 * r["rate"] for r in t]); lo = np.array([100 * r["ci"][0] for r in t])
            hi = np.array([100 * r["ci"][1] for r in t])
            a2.fill_between(xs, lo, hi, color=COLORS[model], alpha=0.12, lw=0)
            a2.plot(xs, y, color=COLORS[model], lw=2, marker="os"[m], ms=6, mfc="white", mew=1.5,
                    label=LABELS[model].replace(" (Pathfinder)", ""))
        a1.axhline(0, color=MUTED, lw=1, ls="--")
        ns = sorted({r["n"] for r in table})
        a1.set_xticks(x, [str(cm) for cm in lengths]); a1.set_xlabel("Genome length (cM)")
        a1.set_ylabel("ΔELBO = best-path ELBO(T_true) − ELBO(T_null) (nats)")
        a1.set_title(f"{SOURCE_LABELS[source]}: how much the true graph is preferred "
                     f"(n = {'–'.join(map(str, ns)) or 0} per box)", loc="left", fontsize=10)
        a1.legend(handles=[Line2D([], [], color=COLORS[m], lw=6, alpha=0.6, label=LABELS[m].replace(" (Pathfinder)", ""))
                           for m in MODEL_NAMES], frameon=False, loc="upper left")
        a2.axhline(50, color=MUTED, lw=1, ls=":"); a2.set_ylim(-3, 103)
        a2.set_xticks(x, [str(cm) for cm in lengths]); a2.set_xlabel("Genome length (cM)")
        a2.set_ylabel("T_true ranked first (%)")
        a2.set_title(f"{SOURCE_LABELS[source]}: correct identification rate (shaded: 95% Wilson interval)",
                     loc="left", fontsize=10)
        a2.legend(frameon=False, loc="lower right")
    fig.suptitle("T_true vs no-admixture tree ((a,b),(c,d)), Pathfinder best-path ELBO", x=0.01, ha="left", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, dest, "topology_rank")


# ---------------------------------------------------------------------------
def parameter_rows(pf, nuts):
    rows = []
    for method, res in (("ibd", pf), ("mixed", pf), ("nuts", nuts)):
        model = "mixed" if method == "nuts" else method
        for (cm, rep, source, m, graph), r in res.items():
            if m != model or graph != "T_true" or r["status"] != "ok":
                continue
            for name, s in r["semantic_parameters"].items():
                rows.append(dict(source=source, method=method, genome_cm=cm, replicate=rep, parameter=name, mean=s["mean"],
                                 lower=s["lower"], upper=s["upper"], truth=s["truth"],
                                 relative_error=(s["mean"] - s["truth"]) / s["truth"]))
    return rows


def parameter_figure(c, rows, dest, source):
    lengths = c["genome_cm"]; names = [p["name"] for p in parameters(c)]
    methods = [m for m in ("ibd", "mixed", "nuts") if any(r["method"] == m for r in rows)]
    if not methods:
        return
    width = 0.78 / len(methods)
    counts = {m: sorted({len({(r["replicate"]) for r in rows if r["method"] == m and r["genome_cm"] == cm})
                         for cm in lengths} - {0}) for m in methods}
    fig, axes = plt.subplots(4, 4, figsize=(19, 13), sharex=True)
    for ax, name in zip(axes.flat, names):
        for k, method in enumerate(methods):
            data = [[r["relative_error"] for r in rows if r["method"] == method and r["parameter"] == name
                     and r["genome_cm"] == cm] for cm in lengths]
            pos = [i for i, d in enumerate(data) if d]
            if not pos:
                continue
            bp = ax.boxplot([data[i] for i in pos], positions=[i + (k - (len(methods) - 1) / 2) * width for i in pos],
                            widths=width * 0.85,
                            patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1.1),
                            whiskerprops=dict(color=MUTED, lw=0.7), capprops=dict(color=MUTED, lw=0.7))
            for box in bp["boxes"]:
                box.set(facecolor=COLORS[method], alpha=0.65, edgecolor=MUTED, lw=0.5)
        truth = next(p["truth"] for p in parameters(c) if p["name"] == name)
        ax.axhline(0, color=INK, lw=1)
        ax.set_title(f"{name}  (truth {truth:,.4g})", fontsize=9.5, loc="left")
        ax.set_xticks(range(len(lengths)), [str(cm) for cm in lengths], fontsize=7)
    for ax in axes[-1]:
        ax.set_xlabel("Genome length (cM)")
    for ax in axes[:, 0]:
        ax.set_ylabel("(estimate − truth) / truth")
    fig.legend(handles=[Line2D([], [], color=COLORS[m], lw=6, alpha=0.65,
                               label=f"{LABELS[m]}, n = {'–'.join(map(str, counts[m])) or 0}") for m in methods],
               loc="upper right", frameon=False, ncol=3, fontsize=9)
    fig.suptitle(f"T_true parameter recovery from {SOURCE_LABELS[source]}: posterior means per subsample (Ne haploid)",
                 x=0.01, ha="left", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    save(fig, dest, f"parameters_{source}")


# ---------------------------------------------------------------------------
def spectrum_figures(c, pf, dest, source):
    cm = c["genome_cm"][-1]; e = edges(c); mid = (e[:-1] + e[1:]) / 2
    pairs = [(i, j) for i in range(len(POPS)) for j in range(i, len(POPS))]
    exposure = pair_counts(n_haploid(c))
    obs, exp, res = defaultdict(list), defaultdict(list), defaultdict(list)
    for (L, rep, _, model, graph), r in pf.items():
        if L != cm or graph != "T_true" or r["status"] != "ok":
            continue
        exp[model].append(np.asarray(r["predictions"]["ibd_number"]) * cm * exposure)
        res[model].append(np.asarray(r["residuals"]["ibd_pearson"]))
        if model == "mixed":
            obs["all"].append(np.asarray(r["observations"]["ibd_count"]))
    if not obs["all"]:
        return
    O = np.array(obs["all"]); n = len(O)
    floor = 0.05          # log axis floor: below one segment, 0 and 0.05 read the same
    fig, axes = plt.subplots(3, 4, figsize=(19, 11), sharex=True)
    for ax, (i, j) in zip(axes.flat, pairs):
        o = O[:, :, i, j]
        med, lo, hi = np.median(o, 0), np.percentile(o, 5, 0), np.percentile(o, 95, 0)
        shown = np.maximum(med, floor)
        ax.errorbar(mid, shown, yerr=[shown - np.maximum(lo, floor), np.maximum(hi, floor) - shown], fmt="none",
                    ecolor=MUTED, elinewidth=1, capsize=0, zorder=3)
        pos = med > 0
        ax.plot(mid[pos], med[pos], "o", ms=4, color=INK, zorder=4, label="observed: median, 5–95%")
        ax.plot(mid[~pos], np.full((~pos).sum(), floor), "o", ms=4, mfc="white", color=INK, zorder=4,
                label="observed median 0 (drawn at floor)")
        for model in MODEL_NAMES:
            if exp[model]:
                E = np.array(exp[model])[:, :, i, j]
                ax.plot(mid, np.median(E, 0), color=COLORS[model], lw=2, label=f"{LABELS[model]}: fitted (median)")
                ax.fill_between(mid, np.percentile(E, 5, 0), np.percentile(E, 95, 0), color=COLORS[model], alpha=0.15, lw=0)
        ax.set_yscale("log"); ax.set_ylim(bottom=floor * 0.8)
        ax.set_title(f"{POPS[i]}–{POPS[j]}", loc="left", fontsize=10)
        ax.grid(axis="y", which="major", color=GRID, lw=0.6); ax.grid(axis="y", which="minor", visible=False)
    for ax in axes.flat[len(pairs):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("Segment length (cM, bin midpoint)")
    for ax in axes[:, 0]:
        ax.set_ylabel("IBD segments per bin")
    seen = {}
    for ax in axes.flat[:len(pairs)]:
        for h, l in zip(*ax.get_legend_handles_labels()):
            seen.setdefault(l, h)
    axes.flat[-1].legend(list(seen.values()), list(seen.keys()), loc="center", frameon=False, fontsize=10)
    fig.suptitle(f"Spectrum fit at {cm} cM, T_true, {SOURCE_LABELS[source]}: every bin, all {n} subsamples "
                 f"(bands: 5–95% across subsamples; axis floored at {floor} segments)",
                 x=0.01, ha="left", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    save(fig, dest, f"spectrum_{cm}cM_{source}")
    fig, axes = plt.subplots(3, 4, figsize=(19, 10), sharex=True, sharey=True)
    for ax, (i, j) in zip(axes.flat, pairs):
        for k, model in enumerate(MODEL_NAMES):
            if res[model]:
                R = np.array(res[model])[:, :, i, j]
                off = (k - 0.5) * 0.12
                ax.errorbar(mid + off, np.median(R, 0), yerr=[np.median(R, 0) - np.percentile(R, 25, 0),
                                                               np.percentile(R, 75, 0) - np.median(R, 0)],
                            fmt="os"[k], ms=4, mfc="white", color=COLORS[model], elinewidth=1.2, capsize=0,
                            label=f"{LABELS[model]}: median, IQR")
        ax.axhline(0, color=INK, lw=1); ax.axhspan(-2, 2, color=GRID, alpha=0.5, lw=0)
        ax.set_title(f"{POPS[i]}–{POPS[j]}", loc="left", fontsize=10)
    for ax in axes.flat[len(pairs):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("Segment length (cM, bin midpoint)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Pearson residual (k − λ)/√λ")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    axes.flat[-1].legend(handles, labels, loc="center", frameon=False, fontsize=10)
    fig.suptitle(f"Residuals at {cm} cM, T_true, {SOURCE_LABELS[source]}, all {n} subsamples (grey band: ±2)",
                 x=0.01, ha="left", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    save(fig, dest, f"residuals_{cm}cM_{source}")


# ---------------------------------------------------------------------------
def hapibd_detection(c, out):
    """hap-IBD over true segment counts per bin and pair, pooled over every block."""
    paths = [out / "pool_000" / f"block_{k:03d}.npz" for k in range(c["blocks"])]
    if not all(p.exists() for p in paths):
        return None
    true = hap = 0
    for p in paths:
        with np.load(p) as z:
            if "hapibd_count" not in z:
                return None
            true = true + z["true_count"]; hap = hap + z["hapibd_count"]
    e = edges(c); rows = []
    for i in range(len(POPS)):
        for j in range(i, len(POPS)):
            for b in range(len(e) - 1):
                t, h = int(true[b, i, j]), int(hap[b, i, j])
                rows.append(dict(bin_lower_cm=e[b], bin_upper_cm=e[b + 1], pair=f"{POPS[i]}-{POPS[j]}", true=t,
                                 hapibd=h, ratio=h / t if t else None))
    iu = np.triu_indices(len(POPS))
    pooled_true = np.array([true[b][iu].sum() for b in range(len(e) - 1)])
    pooled_hap = np.array([hap[b][iu].sum() for b in range(len(e) - 1)])
    return dict(rows=rows, pooled_true=pooled_true, pooled_hap=pooled_hap)


def detection_figure(c, det, dest):
    e = edges(c); mid = (e[:-1] + e[1:]) / 2
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 4.6))
    for pair in sorted({r["pair"] for r in det["rows"]}):
        pts = [(r["bin_lower_cm"], r["ratio"]) for r in det["rows"] if r["pair"] == pair and r["true"] >= 20]
        if pts:
            idx = [list(e).index(x) for x, _ in pts]
            a1.plot(mid[idx], [y for _, y in pts], color=MUTED, lw=0.8, alpha=0.5)
    ok = det["pooled_true"] > 0
    ratio = np.where(ok, det["pooled_hap"] / np.maximum(det["pooled_true"], 1), np.nan)
    a1.plot(mid[ok], ratio[ok], color=INK, lw=2, marker="o", ms=4, label="all pairs pooled")
    a1.plot([], [], color=MUTED, lw=0.8, alpha=0.5, label="one leaf pair (bins with ≥ 20 true segments)")
    a1.axhline(1, color=MUTED, lw=1, ls="--")
    a1.set_xlabel("Segment length (cM, bin midpoint)"); a1.set_ylabel("hap-IBD count / true count")
    a1.set_title("Detection: hap-IBD segments per true segment, by length bin", loc="left", fontsize=10)
    a1.legend(frameon=False, loc="lower right")
    a2.plot(mid, np.maximum(det["pooled_true"], 0.5), color=COLORS["ibd"], lw=2, marker="o", ms=4, label="true IBD")
    a2.plot(mid, np.maximum(det["pooled_hap"], 0.5), color=COLORS["mixed"], lw=2, marker="s", ms=4, mfc="white",
            label="hap-IBD")
    a2.set_yscale("log"); a2.set_xlabel("Segment length (cM, bin midpoint)")
    a2.set_ylabel("segments per bin, all pairs and blocks")
    a2.set_title(f"Pooled spectrum over all {c['blocks']} blocks ({c['blocks'] * c['block_cm']:g} cM)", loc="left",
                 fontsize=10)
    a2.legend(frameon=False, loc="upper right")
    h = c["hapibd"]
    fig.suptitle(f"hap-IBD on the simulated genotypes (min-seed {h['min_seed']:g}, min-output {h['min_output']:g} cM) "
                 "against the true segments", x=0.01, ha="left", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, dest, "hapibd_detection")


def evidence_rows(pf, nuts):
    rows = []
    for (cm, rep, source, model, graph), r in nuts.items():
        p = pf.get((cm, rep, source, model, graph))
        if r["status"] != "ok" or not p or p["status"] != "ok" or "evidence" not in r:
            continue
        rows.append(dict(genome_cm=cm, replicate=rep, source=source, graph=graph, logz=r["evidence"]["logz"],
                         logz_error=r["evidence"]["error"], elbo_max=p["elbo_max"], elbo=p["elbo"],
                         gap=r["evidence"]["logz"] - p["elbo_max"]))
    return rows


def elbo_logz_figure(c, ev, dest):
    lengths = c["genome_cm"]; x = np.arange(len(lengths)); width = 0.36
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(19, 5.2))
    for k, graph in enumerate(GRAPHS):
        data = [[r["gap"] for r in ev if r["graph"] == graph and r["genome_cm"] == cm] for cm in lengths]
        pos = [i for i, d in enumerate(data) if d]
        if pos:
            bp = a1.boxplot([data[i] for i in pos], positions=[i + (k - 0.5) * width for i in pos], widths=width * 0.85,
                            patch_artist=True, showfliers=True, medianprops=dict(color=INK, lw=1.3),
                            flierprops=dict(marker=".", markersize=3, markerfacecolor=MUTED, markeredgecolor=MUTED),
                            whiskerprops=dict(color=MUTED), capprops=dict(color=MUTED))
            for box in bp["boxes"]:
                box.set(facecolor=GRAPH_COLORS[graph], alpha=0.6, edgecolor=MUTED)
    a1.axhline(0, color=MUTED, lw=1, ls="--")
    a1.set_xticks(x, [str(cm) for cm in lengths]); a1.set_xlabel("Genome length (cM)")
    a1.set_ylabel("log Z − best-path ELBO (nats) = KL of the best path")
    a1.set_title("Gap between evidence and ELBO", loc="left", fontsize=10)
    a1.legend(handles=[Line2D([], [], color=GRAPH_COLORS[g], lw=6, alpha=0.6, label=g) for g in GRAPHS],
              frameon=False, loc="upper left")
    paired = defaultdict(dict)
    for r in ev:
        paired[(r["genome_cm"], r["replicate"])][r["graph"]] = r
    pts = [(cm, v["T_true"]["elbo_max"] - v["T_null"]["elbo_max"], v["T_true"]["logz"] - v["T_null"]["logz"],
            np.hypot(v["T_true"]["logz_error"], v["T_null"]["logz_error"]))
           for (cm, rep), v in paired.items() if len(v) == 2]
    if pts:
        cmap = plt.get_cmap("Blues"); cms = np.array([p[0] for p in pts])
        colors = [cmap(0.35 + 0.65 * lengths.index(cm) / max(1, len(lengths) - 1)) for cm in cms]
        de = np.array([p[1] for p in pts]); dz = np.array([p[2] for p in pts])
        a2.scatter(de, dz, c=colors, s=28, edgecolor=INK, lw=0.4, zorder=3)
        lim = [min(de.min(), dz.min(), 0), max(de.max(), dz.max())]
        a2.plot(lim, lim, color=MUTED, lw=1, ls="--", label="ΔlogZ = ΔELBO")
        a2.axhline(0, color=GRID, lw=1); a2.axvline(0, color=GRID, lw=1)
        a2.set_xlabel("ΔELBO(T_true − T_null) (nats)"); a2.set_ylabel("ΔlogZ(T_true − T_null) (nats)")
        a2.set_title("Paired: evidence difference vs ELBO difference", loc="left", fontsize=10)
        a2.legend(handles=[Line2D([], [], color=MUTED, ls="--", label="ΔlogZ = ΔELBO")] +
                  [Line2D([], [], marker="o", ls="", color=cmap(0.35 + 0.65 * i / max(1, len(lengths) - 1)),
                          markeredgecolor=INK, label=f"{cm} cM") for i, cm in enumerate(lengths)],
                  frameon=False, loc="best", fontsize=8, ncol=2)
        # own colours (reference slots 3 and 7): these lines are scores, not graphs
        for label, idx, color, marker in (("by ELBO", 1, "#1baf7a", "o"), ("by log Z", 2, "#4a3aa7", "s")):
            rate, lo, hi, xs = [], [], [], []
            for i, cm in enumerate(lengths):
                g = [p for p in pts if p[0] == cm]
                if g:
                    k = sum(p[idx] > 0 for p in g); w = wilson(k, len(g))
                    xs.append(i); rate.append(100 * k / len(g)); lo.append(100 * w[0]); hi.append(100 * w[1])
            a3.fill_between(xs, lo, hi, color=color, alpha=0.12, lw=0)
            a3.plot(xs, rate, color=color, lw=2, marker=marker, ms=6, mfc="white", mew=1.5, label=label)
        a3.set_ylim(-3, 103); a3.axhline(50, color=MUTED, lw=1, ls=":")
        a3.set_xticks(x, [str(cm) for cm in lengths]); a3.set_xlabel("Genome length (cM)")
        a3.set_ylabel("T_true ranked first (%)")
        a3.set_title(f"Identification by ELBO vs by log Z (n = {len(pts) // max(1, len(set(cms)))} per length)",
                     loc="left", fontsize=10)
        a3.legend(frameon=False, loc="lower right")
    fig.suptitle(f"Mixed model, {SOURCE_LABELS[c['nuts_source']]}: exact evidence (NUTS + bridge sampling) "
                 "vs Pathfinder best-path ELBO", x=0.01,
                 ha="left", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, dest, "elbo_vs_logz")


# ---------------------------------------------------------------------------
def true_graph_figure(c, dest):
    """Merge parents at the midpoint of their children; admixture sources either side of b."""
    d = build(); ne = branch_ne(c)
    level = {"b": 1, "ab": 2, "cb": 3, "cbd": 4, "root": 5}
    x = {"a": 0.4, "b": 1.4, "c": 2.6, "d": 3.6}; y0 = {leaf: 0.0 for leaf in x}; y1 = {}
    for ev in d.ordered_events:
        if ev["type"] == "MERGE":
            for ch in ev["children"]:
                y1[ch] = level[ev["parent"]]
            x[ev["parent"]] = sum(x[ch] for ch in ev["children"]) / 2; y0[ev["parent"]] = level[ev["parent"]]
        else:
            y1[ev["child"]] = level[ev["child"]]
            for s, dx in zip(ev["parents"], (-0.45, 0.45)):
                x[s] = x[ev["child"]] + dx; y0[s] = level[ev["child"]]
    y1["root"] = 5.5
    fig, ax = plt.subplots(figsize=(9, 6.5))
    for node in NODES:
        color = "#eb6834" if node in ("b1", "b2") else "#3d3d3a"
        ax.plot([x[node]] * 2, [y0[node], y1[node]], color=color, lw=1.5 + 5 * ne[node] / max(ne.values()),
                solid_capstyle="butt")
        mid = (y0[node] + min(y1[node], 5.3)) / 2
        ax.text(x[node] + 0.08, mid, f"{node}\nNe {ne[node]:,.0f}", fontsize=8, va="center", bbox=BACK, zorder=6)
    for ev in d.ordered_events:
        if ev["type"] == "MERGE":
            p = ev["parent"]
            for ch in ev["children"]:
                ax.plot([x[ch], x[p]], [y0[p]] * 2, color="#3d3d3a", lw=2)
            ax.scatter([x[p]], [y0[p]], s=26, color=INK, zorder=5)
            ax.text(x[p], y0[p] - 0.1, f"t_{p} = {c['times'][p]:g}", fontsize=8.5, weight="bold", ha="center",
                    va="top", bbox=BACK, zorder=7)
        else:
            ch = ev["child"]
            for s in ev["parents"]:
                ax.plot([x[ch], x[s]], [y1[ch]] * 2, color="#eb6834", lw=2, ls=(0, (3, 1.5)))
            ax.scatter([x[ch]], [y1[ch]], s=26, color=INK, zorder=5)
            ax.text(x[ch], y1[ch] - 0.1, f"t_b_admix = {c['times']['b_admix']:g}", fontsize=8.5, weight="bold",
                    ha="center", va="top", bbox=BACK, zorder=7)
            for s, text in zip(ev["parents"], (f"f_b1 = {c['fraction_b1']:g}", f"1 − f_b1 = {1 - c['fraction_b1']:g}")):
                ax.text(x[s] + 0.08, y1[ch] + 0.22, text, fontsize=8, color="#eb6834", va="center", bbox=BACK, zorder=6)
    for leaf in POPS:
        ax.text(x[leaf], -0.25, leaf, fontsize=14, weight="bold", ha="center", va="top")
    ax.set_yticks(range(6), ["0"] + [f"{c['times'][e]:g}" for e in EVENTS])
    ax.set_ylabel("generations before present (schematic spacing)")
    ax.set_xticks([]); ax.set_xlim(-0.4, 4.4); ax.set_ylim(-0.6, 5.6); ax.grid(False)
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    if "growth" in c:
        g = c["growth"]
        how = (f"growth from {g['ancestral_ne']:,} at the root to present sizes "
               f"{', '.join(f'{k} {v:,}' for k, v in g['present_ne'].items())}")
    else:
        how = "explicit per-branch table"
    ax.set_title(f"Generating graph: haploid Ne per branch ({how}); line width ∝ Ne", fontsize=8.5, loc="left")
    fig.tight_layout()
    save(fig, dest, "topology_true")


# ---------------------------------------------------------------------------
def report(c, dest, topo, params, ev, pf, nuts, detection):
    ne = branch_ne(c)
    lines = ["# Four leaves, one admixture: topology, recovery, ELBO vs log Z", "",
             "Branch Ne (haploid): " + ", ".join(f"{k} {v:,.0f}" for k, v in ne.items()) + ".  IBD sources: "
             + ", ".join(SOURCE_LABELS[s] for s in c["sources"]) + f"; {len(c['bin_edges_cm']) - 1} bins from "
             f"{c['bin_edges_cm'][0]} cM.  Models: mixed and IBD-only Nsmooth (Poisson IBD).  Graph score: "
             f"best-path ELBO.  NUTS on {SOURCE_LABELS[c['nuts_source']]}.", "",
             f"Pathfinder: {len(pf)} fits, {sum(r['status'] != 'ok' for r in pf.values())} failed.  "
             f"NUTS: {len(nuts)} fits, {sum(bool(r.get('reliable')) for r in nuts.values())} pass R-hat/ESS/divergence/treedepth.", "",
             "## T_true vs T_null (Pathfinder)", "",
             "| source | model | cM | n | T_true first | 95% CI | median ΔELBO | smallest ΔELBO |",
             "|---|---|---:|---:|---:|---|---:|---:|"]
    for t in topo["table"]:
        lines.append(f"| {SOURCE_LABELS[t['source']]} | {t['model']} | {t['genome_cm']} | {t['n']} | {t['rate']:.0%} | "
                     f"{t['ci'][0]:.0%}–{t['ci'][1]:.0%} | {t['median']:+.1f} | {t['minimum']:+.1f} |")
    if detection:
        e = edges(c); t, h = detection["pooled_true"], detection["pooled_hap"]
        lines += ["", "## hap-IBD detection (all blocks, all pairs pooled)", "",
                  "| bin (cM) | true | hap-IBD | ratio |", "|---|---:|---:|---:|"]
        lines += [f"| {e[b]:g}–{e[b + 1]:g} | {t[b]} | {h[b]} | " + (f"{h[b] / t[b]:.2f}" if t[b] else "–") + " |"
                  for b in range(len(t))]
    if ev:
        lines += ["", f"## Exact evidence vs best-path ELBO (mixed, {SOURCE_LABELS[c['nuts_source']]}, NUTS subsamples)", "",
                  "| cM | graph | n | median log Z − ELBO | max bridge error |", "|---:|---|---:|---:|---:|"]
        for cm in c["genome_cm"]:
            for g in GRAPHS:
                r = [x for x in ev if x["genome_cm"] == cm and x["graph"] == g]
                if r:
                    lines.append(f"| {cm} | {g} | {len(r)} | {np.median([x['gap'] for x in r]):+.2f} | "
                                 f"{max(x['logz_error'] for x in r):.3f} |")
    columns = [(s, m) for s in c["sources"] for m in ("ibd", "mixed", "nuts") if m != "nuts" or s == c["nuts_source"]]
    names = {"ibd": "IBD-only (PF)", "mixed": "Mixed (PF)", "nuts": "Mixed (NUTS)"}
    lines += ["", "## Parameter recovery at the longest length: median relative error", "",
              "| parameter | truth | " + " | ".join(f"{names[m]}, {SOURCE_LABELS[s]}" for s, m in columns) + " |",
              "|---|---:|" + "---:|" * len(columns)]
    L = c["genome_cm"][-1]
    for p in parameters(c):
        cells = []
        for s, m in columns:
            e = [r["relative_error"] for r in params if r["source"] == s and r["method"] == m and r["genome_cm"] == L
                 and r["parameter"] == p["name"]]
            cells.append(f"{np.median(e):+.0%} (n={len(e)})" if e else "–")
        lines.append(f"| {p['name']} | {p['truth']:,.4g} | " + " | ".join(cells) + " |")
    lines += ["", "Figures: topology_rank, parameters_<source>, spectrum_<L>cM_<source>, residuals_<L>cM_<source>, "
              "hapibd_detection, elbo_vs_logz, topology_true.", ""]
    (dest / "report.md").write_text("\n".join(lines))
