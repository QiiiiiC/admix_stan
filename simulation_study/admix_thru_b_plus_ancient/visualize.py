"""Tables, report and figures.

  contrasts     per subsample: dELBO(T_true - T_alt1) (recent loop kept?) and
                dELBO(T_true - T_alt2) (ancient admixture kept?) for every model
  win_rate      share of subsamples in which T_true beats each alternative, and
                beats BOTH -- the "mixed beats both" claim, as a function of length
  delbo         distribution of each contrast across subsamples, per model
  parameters    relative error of every T_true parameter, per model and length
  topology      the three candidate graphs, events named as in `parameters`

SNP-only fits never read IBD, so the same SNP-only fits appear in both the
true-IBD and the hap-IBD figures.
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

from study import CANDIDATES, CONTRASTS, SOURCES, MODEL_NAMES, save_json

MODELS = list(MODEL_NAMES)
LABELS = {"ibd": "IBD-only", "snp": "SNP-only", "mixed": "Mixed"}
# Reference categorical slots 1-3 (validated all-pairs); identity never by colour alone.
COLORS = {"ibd": "#2a78d6", "snp": "#1baf7a", "mixed": "#eb6834"}
MARKERS = {"ibd": "o", "snp": "s", "mixed": "D"}   # second encoding: coincident lines stay distinguishable
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
CONTRAST_TITLES = {"recent": "T_true − T_alt1: keeps the recent loop?",
                   "ancient": "T_true − T_alt2: keeps the ancient admixture?"}
PARAMETERS = ["loop_split", "loop_donor_merge", "loop_close", "loop_fraction", "ancient_admix",
              "ancient_fraction", "left", "right", "root", "Ne"]

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "axes.grid.axis": "y",
                     "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True})


def csv_write(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def model_source(model, source):
    return "none" if model == "snp" else source


def build_summary(c, out, replicate_limit=None):
    out = Path(out); dest = out / "summary"; (dest / "figures").mkdir(parents=True, exist_ok=True)
    reps = min(c["replicates"], replicate_limit or c["replicates"])
    fits, fit_rows, param_rows = {}, [], []
    base = out / "pool_000" / "fits"
    for path in sorted(base.glob("cm_*/rep_*/*/*/*/result.json")):
        r = json.loads(path.read_text())
        if r["replicate"] >= reps:
            continue
        key = (r["genome_cm"], r["replicate"], r["source"], r["model"], r["candidate"]["id"])
        fits[key] = r
        row = dict(genome_cm=key[0], replicate=key[1], source=key[2], model=key[3], candidate=key[4],
                   status=r["status"], reliable=r["reliable"], elbo=r.get("elbo"), logz=r.get("logz"),
                   ess=r.get("ess"), pareto_k=r.get("pareto_k"), seconds=r["seconds"])
        fit_rows.append(row)
        for name, stat in r.get("semantic_parameters", {}).items():
            param_rows.append(dict(row, parameter=name, mean=stat["mean"], lower=stat["lower"],
                                   upper=stat["upper"], truth=stat["truth"],
                                   relative_error=(stat["mean"] - stat["truth"]) / stat["truth"]))
    contrast_rows = []
    for cm in c["genome_cm"]:
        for rep in range(reps):
            for model in MODELS:
                for source in sorted({model_source(model, s) for s in SOURCES}):
                    got = {k: fits.get((cm, rep, source, model, k)) for k in CANDIDATES}
                    if any(v is None or v["status"] != "ok" for v in got.values()):
                        continue
                    row = dict(genome_cm=cm, replicate=rep, source=source, model=model)
                    for name, (a, b) in CONTRASTS.items():
                        row[f"delbo_{name}"] = got[a]["elbo"] - got[b]["elbo"]
                        row[f"dlogz_{name}"] = got[a]["logz"] - got[b]["logz"]
                    row["true_beats_both"] = row["delbo_recent"] > 0 and row["delbo_ancient"] > 0
                    row["elbo_winner"] = max(CANDIDATES, key=lambda k: got[k]["elbo"])
                    contrast_rows.append(row)
    ident = ["genome_cm", "replicate", "source", "model"]
    csv_write(dest / "fits.csv", fit_rows, ident + ["candidate", "status", "reliable", "elbo", "logz", "ess",
                                                   "pareto_k", "seconds"])
    csv_write(dest / "contrasts.csv", contrast_rows, ident + ["delbo_recent", "delbo_ancient", "dlogz_recent",
                                                             "dlogz_ancient", "true_beats_both", "elbo_winner"])
    csv_write(dest / "parameter_estimates.csv", param_rows, ident + ["candidate", "parameter", "mean", "lower",
                                                                    "upper", "truth", "relative_error"])
    expected = len(c["genome_cm"]) * reps * len(CANDIDATES) * (2 * len(SOURCES) + 1)
    save_json(dest / "summary_status.json", dict(
        subsamples_per_length=reps, expected_fits=expected, attempted=len(fit_rows),
        failed=sum(r["status"] != "ok" for r in fit_rows), reliable=sum(bool(r["reliable"]) for r in fit_rows),
        complete_contrasts=len(contrast_rows),
        note="Subsamples overlap (blocks drawn without replacement from one pool of "
             f"{c['blocks']}); ELBO is the ranking statistic."))
    table = rates(c, contrast_rows)
    report(c, dest, table, fit_rows, reps)
    for source in SOURCES:
        rows = [r for r in contrast_rows if r["source"] == model_source(r["model"], source)]
        win_rate_figure(c, [t for t in table if t["source"] == model_source(t["model"], source)], source, dest)
        delbo_figure(c, rows, source, dest)
        prm = [r for r in param_rows if r["candidate"] == "T_true" and r["source"] == model_source(r["model"], source)]
        parameter_figure(c, prm, source, dest)
    from plot_topologies import topology_figure
    topology_figure(c, dest / "figures")
    print(f"Summarised {len(fit_rows)} fits ({reps} subsamples per length) into {dest}")


def rates(c, contrast_rows):
    groups = defaultdict(list)
    for r in contrast_rows:
        groups[(r["source"], r["model"], r["genome_cm"])].append(r)
    table = []
    for (source, model, cm), rows in sorted(groups.items()):
        row = dict(source=source, model=model, genome_cm=cm, n=len(rows),
                   both=float(np.mean([r["true_beats_both"] for r in rows])))
        for name in CONTRASTS:
            d = np.array([r[f"delbo_{name}"] for r in rows])
            row.update({f"median_{name}": float(np.median(d)), f"win_{name}": float(np.mean(d > 0))})
        table.append(row)
    return table


def report(c, dest, table, fit_rows, reps):
    t = c["times"]; f = c["fractions"]
    lines = ["# Recent loop + ancient admixture (Nfixed): IBD-only vs SNP-only vs Mixed", "",
             f"Truth: loop on a {t['loop_split']}→{t['loop_donor_merge']}→{t['loop_close']} gen "
             f"(α = {f['loop']}), ancient admixture on c at {t['ancient_admix']} gen ({f['ancient']:.0%}/"
             f"{1-f['ancient']:.0%}), merges {t['left']}/{t['right']}/{t['root']}, haploid Ne "
             f"{c['haploid_ne']:,.0f} everywhere.  {reps} subsamples per genome length, blocks of "
             f"{c['block_cm']:g} cM drawn without replacement from a pool of {c['blocks']}; subsamples "
             "overlap, so shares are sensitivities to genome selection, not independent-replicate rates.", "",
             "ΔELBO = ELBO(T_true) − ELBO(alternative), each from three pooled Pathfinder starts on normalised "
             "models. 'wins' = share of subsamples with ΔELBO > 0; 'both' = share in which T_true beats both "
             "alternatives. SNP-only is fitted once per subsample and repeated under both IBD sources.", ""]
    n_ok = sum(r["status"] == "ok" for r in fit_rows)
    lines += [f"Fits: {len(fit_rows)} attempted, {n_ok} ok, "
              f"{sum(bool(r['reliable']) for r in fit_rows)} pass the ESS/Pareto-k screen.", ""]
    for source in SOURCES:
        lines += [f"## IBD source: {source}", "",
                  "| cM | model | median Δ recent | wins recent | median Δ ancient | wins ancient | both | n |",
                  "|---:|---|---:|---:|---:|---:|---:|---:|"]
        for cm in c["genome_cm"]:
            for model in MODELS:
                row = next((x for x in table if x["genome_cm"] == cm and x["model"] == model
                            and x["source"] == model_source(model, source)), None)
                if row:
                    lines.append(f"| {cm} | {LABELS[model]} | {row['median_recent']:+.1f} | {row['win_recent']:.0%} | "
                                 f"{row['median_ancient']:+.1f} | {row['win_ancient']:.0%} | {row['both']:.0%} | {row['n']} |")
        lines.append("")
    lines += ["Figures: `figures/win_rate_*`, `figures/delbo_*`, `figures/parameters_*` (one per IBD source).",
              "Tables: `contrasts.csv` (per subsample), `fits.csv`, `parameter_estimates.csv`.", ""]
    (dest / "report.md").write_text("\n".join(lines))


def save(fig, dest, name):
    for ext in ("png", "pdf"):
        fig.savefig(dest / "figures" / f"{name}.{ext}", dpi=200 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def length_axis(ax, lengths):
    ax.set_xticks(range(len(lengths)), [str(x) for x in lengths])
    ax.set_xlabel("Genome length (cM)")


def spread(values, gap):
    """Label positions no closer than `gap`, as near the targets as possible (sorted sweep, then recentred)."""
    order = np.argsort(values)
    placed = np.asarray(values, float)[order].copy()
    for k in range(1, len(placed)):
        placed[k] = max(placed[k], placed[k - 1] + gap)
    placed -= (placed - np.asarray(values, float)[order]).mean()
    out = np.empty_like(placed); out[order] = placed
    return out


def win_rate_figure(c, table, source, dest):
    lengths = c["genome_cm"]
    panels = [("win_recent", CONTRAST_TITLES["recent"]), ("win_ancient", CONTRAST_TITLES["ancient"]),
              ("both", "T_true beats BOTH alternatives")]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    for ax, (key, title) in zip(axes, panels):
        ends = []
        for model in MODELS:
            rows = {r["genome_cm"]: r for r in table if r["model"] == model}
            x = [i for i, cm in enumerate(lengths) if cm in rows]
            y = [100 * rows[lengths[i]][key] for i in x]
            if not x:
                continue
            ax.plot(x, y, color=COLORS[model], lw=2, marker=MARKERS[model], ms=6, mfc="white", mew=1.5,
                    label=LABELS[model])
            ends.append((model, x[-1], y[-1]))
        for (model, x, y), ly in zip(ends, spread([e[2] for e in ends], 7.0)):
            ax.annotate(LABELS[model], (x, y), xytext=(x + 0.25, ly), textcoords="data",
                        va="center", fontsize=8, color=INK,
                        arrowprops=dict(arrowstyle="-", color=COLORS[model], lw=0.8) if abs(ly - y) > 1 else None)
        ax.set_title(title, fontsize=9.5, color=INK, loc="left")
        ax.set_ylim(-8, 108); ax.set_xlim(-0.3, len(lengths) - 0.3 + 1.2)
        length_axis(ax, lengths)
    axes[0].set_ylabel("Subsamples preferring T_true (%)")
    axes[0].legend(frameon=False, loc="lower right", fontsize=8)
    fig.suptitle(f"How often each model prefers the true graph — {source} IBD", x=0.01, ha="left",
                 fontsize=11, color=INK)
    fig.tight_layout()
    save(fig, dest, f"win_rate_{source}")


def delbo_figure(c, rows, source, dest):
    lengths = c["genome_cm"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 7.5))
    for i, name in enumerate(CONTRASTS):
        for j, model in enumerate(MODELS):
            ax = axes[i, j]
            data = [[r[f"delbo_{name}"] for r in rows if r["model"] == model and r["genome_cm"] == cm]
                    for cm in lengths]
            pos = [k for k, d in enumerate(data) if d]
            if pos:
                bp = ax.boxplot([data[k] for k in pos], positions=pos, widths=0.55, patch_artist=True,
                                showfliers=False, medianprops=dict(color=INK, lw=1.5),
                                whiskerprops=dict(color=MUTED), capprops=dict(color=MUTED))
                for box in bp["boxes"]:
                    box.set(facecolor=COLORS[model], alpha=0.55, edgecolor=MUTED)
            ax.axhline(0, color=MUTED, lw=1, ls="--")
            ax.set_xlim(-0.6, len(lengths) - 0.4)
            ax.set_xticks(range(len(lengths)), [f"{cm}\n{np.mean(np.array(d) > 0):.0%}" if d else str(cm)
                                                for cm, d in zip(lengths, data)], fontsize=8)
            ax.set_xlabel("Genome length (cM) / share of subsamples with ΔELBO > 0")
            ax.set_title(f"{LABELS[model]}: {CONTRAST_TITLES[name]}", fontsize=9, color=INK, loc="left")
            if j == 0:
                ax.set_ylabel("ΔELBO (nats)")
    fig.suptitle(f"ΔELBO across subsamples, {source} IBD (each panel has its own scale)",
                 x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout()
    save(fig, dest, f"delbo_{source}")


def subsample_counts(rows, lengths):
    """Number of subsamples behind each box.  The title always states it; tick
    labels carry it (IBD/SNP/Mixed in order where models differ) only when it is
    not the same for every box, so uniform counts do not crowd the axis."""
    per = {cm: [len({r["replicate"] for r in rows if r["model"] == m and r["genome_cm"] == cm})
                for m in MODELS] for cm in lengths}
    counts = [n for v in per.values() for n in v]
    uniform = len(set(counts)) <= 1
    labels = [str(cm) if uniform else
              f"{cm}\nn={v[0]}" if len(set(v)) == 1 else f"{cm}\nn={'/'.join(map(str, v))}"
              for cm, v in per.items()]
    return labels, uniform, (min(counts), max(counts)) if counts else (0, 0)


def parameter_figure(c, rows, source, dest):
    lengths = c["genome_cm"]
    ticks, uniform, (n_lo, n_hi) = subsample_counts(rows, lengths)
    fig, axes = plt.subplots(2, 5, figsize=(18, 6.8), sharex=True)
    width = 0.26
    for ax, name in zip(axes.flat, PARAMETERS):
        for m, model in enumerate(MODELS):
            data = [[r["relative_error"] for r in rows if r["model"] == model and r["parameter"] == name
                     and r["genome_cm"] == cm] for cm in lengths]
            pos = [k for k, d in enumerate(data) if d]
            if not pos:
                continue
            bp = ax.boxplot([data[k] for k in pos], positions=[k + (m - 1) * width for k in pos], widths=width * 0.85,
                            patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1.2),
                            whiskerprops=dict(color=MUTED, lw=0.8), capprops=dict(color=MUTED, lw=0.8))
            for box in bp["boxes"]:
                box.set(facecolor=COLORS[model], alpha=0.6, edgecolor=MUTED, lw=0.6)
        ax.axhline(0, color=INK, lw=1)
        ax.set_title(name, fontsize=9.5, color=INK, loc="left")
        ax.set_xticks(range(len(lengths)), ticks, fontsize=7)
    for ax in axes[1]:
        ax.set_xlabel("Genome length (cM)" if uniform else "Genome length (cM) / subsamples per box")
    for ax in axes[:, 0]:
        ax.set_ylabel("(estimate − truth) / truth")
    handles = [Line2D([], [], color=COLORS[m], lw=6, alpha=0.6, label=LABELS[m]) for m in MODELS]
    fig.legend(handles=handles, loc="upper right", frameon=False, ncol=3, fontsize=9)
    n_text = f"{n_lo}" if n_lo == n_hi else f"{n_lo}–{n_hi}"
    fig.suptitle(f"T_true parameter recovery, {source} IBD: posterior means, n = {n_text} subsamples per box "
                 "(Ne = the single shared effective_N of the Nfixed models; SNP-only has no IBD, same in both)",
                 x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, dest, f"parameters_{source}")
