"""Tables and figures.  Beyond topology recovery and the spectrum fit:

  evidence_benefit     paired logZ / ELBO (separate - shared) on the SAME data
                       and the SAME graph -- the direct test of whether a second
                       Ne trajectory is supported once its extra parameters are
                       paid for (constants are normalised in run_study).
  ne_trajectory        the true effective Ne(t) of every backbone branch with
                       the fitted constants laid over it, plus the harmonic mean
                       (what accumulated drift implies) and the arithmetic mean.
  ne_ibd_vs_snp        log(Ne_IBD / Ne_SNP) per branch from the separate model.
                       Growth is expected to push this ABOVE zero on the growing
                       branches: IBD sees the recent, larger end; SNP drift sees
                       the harmonic mean, which is dominated by the small old end.
  residual_tilt        slope of the Pearson residual against segment length --
                       the fingerprint a constant Ne leaves when the truth grew.
  component_likelihood lp and chi2/n by data type, shared vs separate.

Ne and event times are NOT expected to be recovered under `growth`: the truth
is outside every candidate's model class.  Those plots are descriptive.
"""
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from scipy.stats import chi2 as chi2_dist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from study import (SCENARIOS, SOURCES, VARIANTS, BACKBONE, candidates, edges, pair_counts,
                   reference, effective_ne, branch_map, save_json)

PAIRS = [(i, j) for i in range(3) for j in range(i, 3)]
PAIR_NAMES = ["a–a", "a–b", "a–c", "b–b", "b–c", "c–c"]
COLORS = ["#236FA6", "#D87921", "#269B79", "#925BB7"]
LABELS = ["Omit\nshared", "Omit\nseparate", "Explicit\nshared", "Explicit\nseparate"]
COMPARISONS = [("omitted_backbone", "poisson_shared"), ("omitted_backbone", "poisson_separate"),
               ("explicit_loop", "poisson_shared"), ("explicit_loop", "poisson_separate")]
NE_COMPONENTS = [("poisson_shared", "Ne"), ("poisson_separate", "Ne_ibd"), ("poisson_separate", "Ne_snp")]
NE_LABELS = ["Shared", "Separate\nIBD", "Separate\nSNP"]
NE_COLORS = ["#236FA6", "#D87921", "#269B79"]


def csv_write(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def poisson_band(k, level=0.6827):
    al = 1-level; k = np.asarray(k, float)
    lo = np.where(k > 0, 0.5*chi2_dist.ppf(al/2, np.maximum(2*k, 1e-9)), 0.0)
    hi = 0.5*chi2_dist.ppf(1-al/2, 2*k+2)
    return lo, hi


def build_summary(c, out):
    out = Path(out); dest = out / "summary"; dest.mkdir(parents=True, exist_ok=True)
    graphs = defaultdict(list)
    for cand in candidates(): graphs[cand["graph"]].append(cand)
    ranking_rows = []; selection_rows = []; diagnostics = []; parameters = []; ne_rows = []
    tilt_rows = []; lik_rows = []; index_rows = []
    plotted = defaultdict(list); score_table = {}
    for scenario in SCENARIOS:
        for pool in range(c["pools"]):
            for rep in range(c["replicates"]):
                for source in SOURCES:
                    for variant in VARIANTS:
                        ident = dict(scenario=scenario, pool=pool, replicate=rep, source=source, variant=variant)
                        folder = out / scenario / f"pool_{pool:03d}" / "fits" / f"rep_{rep:03d}" / source / variant
                        fits = {p.parent.name: json.loads(p.read_text()) for p in folder.glob("*/result.json")}
                        scores = {}; best = {}; reliability = {}
                        for gid, orders in graphs.items():
                            rows = [fits.get(o["id"]) for o in orders]
                            if any(r is None or r["status"] != "ok" for r in rows): continue
                            scores[gid] = dict(logz=float(logsumexp([r["logz"] for r in rows])-np.log(len(rows))),
                                               elbo=float(logsumexp([r["elbo"] for r in rows])-np.log(len(rows))))
                            best[gid] = max(rows, key=lambda r: r["elbo"])
                            reliability[gid] = all(r["reliable"] for r in rows)
                        score_table[(scenario, pool, rep, source, variant)] = scores
                        for cid, r in fits.items():
                            g = r["candidate"]["graph"]
                            # best_order marks the event order the figures use for this graph on
                            # this subsample; the report's per-branch tables filter on it so they
                            # never pool two orders' sizes into one median.
                            base = dict(ident, candidate=cid, graph=g, reliable=r["reliable"],
                                        best_order=(g in best and best[g]["candidate"]["id"] == cid))
                            diagnostics.append(dict(base, status=r["status"], ess=r.get("ess"),
                                                    pareto_k=r.get("pareto_k"), seconds=r["seconds"]))
                            index_rows.append(dict(base, result=str((folder/cid/"result.json").relative_to(out))))
                            if r["status"] != "ok": continue
                            for name, stat in r.get("semantic_parameters", {}).items():
                                parameters.append(dict(base, parameter=name, mean=stat["mean"], lower=stat["lower"],
                                    upper=stat["upper"], truth=stat["truth"], error=stat["mean"]-stat["truth"]))
                            for key, component in (("Ne", "shared"), ("Ne_ibd", "IBD"), ("Ne_snp", "SNP")):
                                if key not in r: continue
                                for j, (node, ref) in enumerate(zip(r["candidate"]["nodes"], r["ne_reference"])):
                                    ne_rows.append(dict(base, node=node, branch=ref.get("branch"), component=component,
                                        mean=r[key]["mean"][j], lower=r[key]["lower"][j], upper=r[key]["upper"][j],
                                        **{k: ref.get(k) for k in ("t0", "t1", "young", "old", "harmonic", "arithmetic")}))
                            for pair, stat in r["residual_tilt"].items():
                                tilt_rows.append(dict(base, pair=pair, **stat))
                            lik = r.get("likelihood_summaries", {})
                            lik_rows.append(dict(base, **{k: lik[k]["mean"] for k in lik}))
                        for search, ids in (("loop_omitted", list(range(1, 22))), ("including_loop", list(range(1, 23)))):
                            complete = all(g in scores for g in ids)
                            row = dict(ident, search=search, complete=complete,
                                       reliable=complete and all(reliability[g] for g in ids),
                                       winner=None, backbone_correct=None, full_graph_correct=None, elbo_gap=None,
                                       rank_backbone=None, rank_explicit=None, logz_winner=None)
                            if complete:
                                # ELBO is the primary ranking statistic.  The importance-sampling
                                # logZ is kept as a secondary column only: pooled ESS is O(1-20)
                                # of 12,000 draws and Pareto k > 1 on essentially every fit, so
                                # that estimate is one draw's weight, not a marginal likelihood.
                                ordered = sorted(ids, key=lambda g: scores[g]["elbo"], reverse=True)
                                win = ordered[0]
                                row.update(winner=win, backbone_correct=win in (15, 22), full_graph_correct=win == 22,
                                           elbo_gap=scores[win]["elbo"]-scores[ordered[1]]["elbo"],
                                           rank_backbone=ordered.index(15)+1,
                                           rank_explicit=(ordered.index(22)+1 if 22 in ids else None),
                                           logz_winner=max(ids, key=lambda g: scores[g]["logz"]))
                                plotted[(scenario, source, "winner_"+search, variant)].append(best[win])
                                for rank, g in enumerate(ordered, 1):
                                    ranking_rows.append(dict(ident, search=search, rank=rank, graph=g,
                                        newick=graphs[g][0]["newick"], best_order=best[g]["candidate"]["id"],
                                        reliable=reliability[g], **scores[g]))
                            selection_rows.append(row)
                        for g, label in ((15, "omitted_backbone"), (22, "explicit_loop")):
                            if g in best: plotted[(scenario, source, label, variant)].append(best[g])
    # Paired shared-vs-separate on identical data.
    paired = []
    for (scenario, pool, rep, source, variant), sc in score_table.items():
        if variant != "poisson_shared": continue
        other = score_table.get((scenario, pool, rep, source, "poisson_separate"), {})
        row = dict(scenario=scenario, pool=pool, replicate=rep, source=source)
        for g, tag in ((15, "backbone"), (22, "explicit")):
            row[f"dlogz_{tag}"] = other[g]["logz"]-sc[g]["logz"] if g in sc and g in other else None
            row[f"delbo_{tag}"] = other[g]["elbo"]-sc[g]["elbo"] if g in sc and g in other else None
        sel = {r["variant"]: r for r in selection_rows if r["scenario"] == scenario and r["pool"] == pool
               and r["replicate"] == rep and r["source"] == source and r["search"] == "including_loop"}
        if all(v in sel and sel[v]["complete"] for v in VARIANTS):
            row.update(shared_backbone=sel["poisson_shared"]["backbone_correct"],
                       separate_backbone=sel["poisson_separate"]["backbone_correct"],
                       shared_winner=sel["poisson_shared"]["winner"], separate_winner=sel["poisson_separate"]["winner"])
        paired.append(row)
    common = ["scenario", "pool", "replicate", "source", "variant"]
    csv_write(dest/"search_rankings.csv", ranking_rows, common+["search", "rank", "graph", "newick", "best_order", "reliable", "logz", "elbo"])
    csv_write(dest/"search_selection.csv", selection_rows, common+["search", "complete", "reliable", "winner", "backbone_correct", "full_graph_correct", "elbo_gap", "rank_backbone", "rank_explicit", "logz_winner"])
    csv_write(dest/"paired_ne_comparison.csv", paired, ["scenario", "pool", "replicate", "source", "dlogz_backbone", "delbo_backbone", "dlogz_explicit", "delbo_explicit", "shared_backbone", "separate_backbone", "shared_winner", "separate_winner"])
    csv_write(dest/"fit_diagnostics.csv", diagnostics, common+["candidate", "graph", "status", "reliable", "ess", "pareto_k", "seconds"])
    csv_write(dest/"parameter_estimates.csv", parameters, common+["candidate", "graph", "reliable", "parameter", "mean", "lower", "upper", "truth", "error"])
    csv_write(dest/"ne_estimates.csv", ne_rows, common+["candidate", "graph", "reliable", "best_order", "node", "branch", "component", "mean", "lower", "upper", "t0", "t1", "young", "old", "harmonic", "arithmetic"])
    csv_write(dest/"residual_tilt.csv", tilt_rows, common+["candidate", "graph", "reliable", "best_order", "pair", "slope", "correlation"])
    csv_write(dest/"component_likelihood.csv", lik_rows, common+["candidate", "graph", "reliable", "best_order", "lp_ibd", "lp_snp", "chi2_ibd", "chi2_snp"])
    save_json(dest/"fit_index.json", index_rows)
    save_json(dest/"summary_status.json", dict(expected_fits=29*len(VARIANTS)*len(SOURCES)*len(SCENARIOS)*c["replicates"]*c["pools"],
        attempted=len(diagnostics), failed=sum(r["status"] != "ok" for r in diagnostics),
        reliable=sum(r["reliable"] for r in diagnostics),
        note="Overlapping genome subsamples; composite-likelihood importance estimates; descriptive plots include unreliable fits."))
    report(c, dest, selection_rows, paired, ne_rows, tilt_rows)
    figures(c, plotted, paired, selection_rows, dest)
    print(f"Archived {len(diagnostics)} fits; tables and figures in {dest}")


def report(c, dest, selection_rows, paired, ne_rows, tilt_rows):
    L = ["# Growth + asymmetric-Ne b-loop: shared vs separate Ne", "",
         "Overlapping 30-of-50 block subsamples; no independent-replicate error bars. Evidence is compared only",
         "within a model variant and IBD source (graph scores average over event orders). Under `growth` every",
         "candidate is misspecified in Ne, so times and sizes are descriptive, not recovery targets.", "",
         "## Topology selection (ELBO ranking, search including the explicit-loop graph)", "",
         "Ranked by ELBO averaged over event orders. The importance-sampling logZ is recorded but not used to rank:",
         "pooled ESS is O(1–20) of 12,000 draws and Pareto k > 1 on essentially every fit, so it is one draw's weight.", "",
         "| Scenario | IBD | Model | Complete | Reliable | Backbone correct | Explicit loop wins | Median rank of g15 | Median rank of g22 |",
         "|---|---|---|---:|---:|---:|---:|---:|---:|"]
    for scenario in SCENARIOS:
        for source in SOURCES:
            for variant in VARIANTS:
                rows = [r for r in selection_rows if r["scenario"] == scenario and r["source"] == source
                        and r["variant"] == variant and r["search"] == "including_loop" and r["complete"]]
                n = len(rows)
                if not n:
                    L.append(f"| {scenario} | {source} | {variant} | 0 | 0 | — | — | — | — |"); continue
                L.append(f"| {scenario} | {source} | {variant} | {n} | {sum(r['reliable'] for r in rows)} | "
                         f"{np.mean([r['backbone_correct'] for r in rows]):.0%} | {np.mean([r['full_graph_correct'] for r in rows]):.0%} | "
                         f"{np.median([r['rank_backbone'] for r in rows]):.0f} | {np.median([r['rank_explicit'] for r in rows]):.0f} |")
    L += ["", "## Does a separate Ne help?  Paired on identical data", "",
          "dELBO = separate − shared for the same graph on the same subsample; positive favours separate. Constants",
          "are normalised so the extra parameters are paid for. 'Rescue' = separate picks a correct backbone where",
          "shared did not; 'loss' the reverse. (dlogZ is in paired_ne_comparison.csv, subject to the caveat above.)", "",
          "| Scenario | IBD | n | mean dELBO backbone | share > 0 | mean dELBO explicit | share > 0 | rescue | loss |",
          "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for scenario in SCENARIOS:
        for source in SOURCES:
            rows = [r for r in paired if r["scenario"] == scenario and r["source"] == source and r["delbo_backbone"] is not None]
            if not rows: continue
            b = np.array([r["delbo_backbone"] for r in rows]); e = np.array([r["delbo_explicit"] for r in rows if r["delbo_explicit"] is not None])
            sel = [r for r in rows if "shared_backbone" in r]
            rescue = sum((not r["shared_backbone"]) and r["separate_backbone"] for r in sel)
            loss = sum(r["shared_backbone"] and not r["separate_backbone"] for r in sel)
            L.append(f"| {scenario} | {source} | {len(rows)} | {b.mean():+.1f} | {np.mean(b > 0):.0%} | "
                     f"{(e.mean() if len(e) else float('nan')):+.1f} | {(np.mean(e > 0) if len(e) else float('nan')):.0%} | {rescue} | {loss} |")
    L += ["", "## Branch sizes on the omitted-loop backbone (graph 15), median fitted vs truth", "",
          "Truth columns are the true effective trajectory: size at the young end, old end, its harmonic mean",
          "(accumulated drift ⇔ duration / harmonic) and arithmetic mean. For b the truth is the composite of the",
          "leaf, the loop (1 / Σ p_k²/N_k) and anc_b. Fitted = median of posterior means across subsamples, at the", "best-ELBO event order of g15 on each subsample (the same fits the figures use).", ""]
    for scenario in SCENARIOS:
        L += [f"### {scenario}", "", "| branch | young | old | harmonic | arithmetic | " +
              " | ".join(f"{s} {lab}" for s in SOURCES for lab in ("shared", "sep IBD", "sep SNP")) + " |",
              "|---|---:|---:|---:|---:|" + "---:|"*6]
        for branch in BACKBONE:
            ref = reference(c, scenario, branch, collapsed=True)
            cells = []
            for source in SOURCES:
                for variant, comp in (("poisson_shared", "shared"), ("poisson_separate", "IBD"), ("poisson_separate", "SNP")):
                    v = [r["mean"] for r in ne_rows if r["scenario"] == scenario and r["source"] == source and r["variant"] == variant
                         and r["graph"] == 15 and r["best_order"] and r["component"] == comp and r["branch"] == branch]
                    cells.append(f"{np.median(v):,.0f}" if v else "—")
            L.append(f"| {branch} | {ref['young']:,.0f} | {ref['old']:,.0f} | {ref['harmonic']:,.0f} | {ref['arithmetic']:,.0f} | " + " | ".join(cells) + " |")
        L.append("")
    L += ["## Residual tilt on the backbone (Pearson residual per cM), median across subsamples", "",
          "| Scenario | IBD | Model | " + " | ".join(PAIR_NAMES) + " |", "|---|---|---|" + "---:|"*6]
    for scenario in SCENARIOS:
        for source in SOURCES:
            for variant in VARIANTS:
                cells = []
                for i, j in PAIRS:
                    v = [r["slope"] for r in tilt_rows if r["scenario"] == scenario and r["source"] == source
                         and r["variant"] == variant and r["graph"] == 15 and r["best_order"] and r["pair"] == f"{i},{j}"]
                    cells.append(f"{np.median(v):+.3f}" if v else "—")
                L.append(f"| {scenario} | {source} | {variant} | " + " | ".join(cells) + " |")
    L += ["", "Figures: `figures/` — topology_recovery, evidence_benefit, spectrum_fit_*, ne_trajectory_*, ne_ibd_vs_snp_*,",
          "residual_tilt_*, component_likelihood_*, parameters_backbone_*, topology_generating_*."]
    (dest/"report.md").write_text("\n".join(L)+"\n")


def box(ax, groups, truth=None, labels=LABELS, colors=COLORS):
    for j, group in enumerate(groups):
        if not group: continue
        bp = ax.boxplot([group], positions=[j], widths=.55, patch_artist=True, showfliers=False,
                        medianprops=dict(color="#222222"))
        bp["boxes"][0].set(facecolor=colors[j], alpha=.35)
        ax.scatter(j+np.linspace(-.14, .14, len(group)), group, s=14, color=colors[j], zorder=3)
    if truth is not None: ax.axhline(truth, color="#A42E32", ls="--", lw=1.1)
    ax.set_xticks(range(len(groups)), [f"{labels[i]}\nn={len(g)}" for i, g in enumerate(groups)], fontsize=7)
    ax.set_xlim(-.6, len(groups)-.4); ax.grid(axis="y", alpha=.15)
    ax.spines[["top", "right"]].set_visible(False)


def save(fig, folder, name, subtitle):
    fig.text(.5, .015, subtitle, ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .06, 1, .95))
    fig.savefig(folder/(name+".png"), dpi=160); fig.savefig(folder/(name+".pdf")); plt.close(fig)


def figures(c, data, paired, selection_rows, dest):
    folder = dest/"figures"; folder.mkdir(exist_ok=True)
    from plot_topologies import topology_figures
    topology_figures(c, folder)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})

    # ---- 1. topology recovery --------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
    for ax, scenario in zip(axes, SCENARIOS):
        k = 0
        for source in SOURCES:
            for variant in VARIANTS:
                rows = [r for r in selection_rows if r["scenario"] == scenario and r["source"] == source
                        and r["variant"] == variant and r["search"] == "including_loop" and r["complete"]]
                if rows:
                    ax.bar(k-.2, np.mean([r["backbone_correct"] for r in rows]), .38, color="#236FA6")
                    ax.bar(k+.2, np.mean([r["full_graph_correct"] for r in rows]), .38, color="#D87921")
                    ax.text(k, 1.02, f"n={len(rows)}", ha="center", fontsize=8)
                k += 1
        ax.set_xticks(range(4), [f"{s}\n{v.split('_')[1]}" for s in SOURCES for v in VARIANTS], fontsize=8)
        ax.set_ylim(0, 1.12); ax.set_title(scenario); ax.grid(axis="y", alpha=.15)
    axes[0].set_ylabel("frequency across subsamples")
    fig.legend([Line2D([0], [0], color="#236FA6", lw=8), Line2D([0], [0], color="#D87921", lw=8)],
               ["backbone ((a,b.1),(b.2,c)) selected (g15 or g22)", "explicit-loop graph g22 selected"],
               loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(.5, .97))
    save(fig, folder, "topology_recovery", "Search over all 22 graphs; ELBO averaged over event orders. Overlapping subsamples.")

    # ---- 2. evidence benefit ---------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5))
    for col, scenario in enumerate(SCENARIOS):
        for row, key in enumerate(("delbo", "dlogz")):
            groups, labels = [], []
            for tag in ("backbone", "explicit"):
                for source in SOURCES:
                    groups.append([r[f"{key}_{tag}"] for r in paired if r["scenario"] == scenario and r["source"] == source
                                   and r[f"{key}_{tag}"] is not None])
                    labels.append(f"g{15 if tag=='backbone' else 22}\n{source}")
            ax = axes[row, col]
            box(ax, groups, 0.0, labels, ["#236FA6", "#236FA6", "#925BB7", "#925BB7"])
            ax.set_title(f"{scenario}: {'logZ' if key=='dlogz' else 'ELBO'} (separate − shared)")
            ax.set_ylabel("nats; > 0 favours separate Ne")
    save(fig, folder, "evidence_benefit",
         "Same subsample, same graph, two models; constants normalised so parameter count is paid for. Points = subsamples.")

    for scenario in SCENARIOS:
        for source in SOURCES:
            # ---- 3. spectrum fit ------------------------------------------------
            x = (edges(c)[:-1]+edges(c)[1:])/2
            ex = next((r for lab, v in COMPARISONS for r in data[(scenario, source, lab, v)] if r["replicate"] == 0), None)
            if ex is not None:
                fig, axes = plt.subplots(2, 6, figsize=(23, 7.4), sharex=True)
                counts = np.asarray(ex["observations"]["ibd_count"]); cm = ex["observations"]["cm"]
                npair = np.asarray(pair_counts(c))
                unit = 1.0
                lo, hi = poisson_band(counts)
                for p, (i, j) in enumerate(PAIRS):
                    k = counts[:, i, j]; m = k > 0
                    axes[0, p].errorbar(x[m], k[m], yerr=[k[m]-lo[:, i, j][m], hi[:, i, j][m]-k[m]], fmt="o", ms=3,
                                        lw=.8, color="k", zorder=3, label="observed (rep 0)")
                    if (~m).any():
                        axes[0, p].plot(x[~m], np.full((~m).sum(), 2.996*unit), "v", ms=3, color=".6", label="0 (95% UL)")
                    for q, (lab, v) in enumerate(COMPARISONS):
                        r0 = next((r for r in data[(scenario, source, lab, v)] if r["replicate"] == 0), None)
                        if r0 is None: continue
                        lam = np.asarray(r0["predictions"]["ibd_number"])[:, i, j]*cm*npair[i, j]
                        axes[0, p].plot(x, lam, color=COLORS[q], lw=1.5, label=LABELS[q].replace("\n", " "))
                    # Floor at a tenth of one segment, ceiling at 10x the largest count in the
                    # whole figure (never 0, which would invert the axis on an empty pair).
                    axes[0, p].set_yscale("log"); axes[0, p].set_ylim(0.1*unit, max(counts.max(), 1)*10)
                    axes[0, p].set_title(f"{PAIR_NAMES[p]}  n={k.sum():,}", fontsize=10)
                    ax = axes[1, p]
                    for q, (lab, v) in enumerate(COMPARISONS):
                        arr = np.asarray([r["residuals"]["ibd_pearson"] for r in data[(scenario, source, lab, v)]])
                        if len(arr):
                            y = arr[:, :, i, j]; qs = np.quantile(y, [.25, .5, .75], axis=0)
                            ax.plot(x, qs[1], color=COLORS[q], lw=1.3); ax.fill_between(x, qs[0], qs[2], color=COLORS[q], alpha=.12)
                    ax.axhline(0, color="k", lw=.8)
                    for yy in (-2, 2): ax.axhline(yy, color="grey", lw=.6, ls=":")
                    ax.set_xlabel("segment length (cM)")
                axes[0, 0].set_ylabel("segments per bin"); axes[1, 0].set_ylabel("Pearson residual (k−λ)/√λ")
                axes[0, 0].legend(fontsize=7)
                fig.suptitle(f"{scenario} / {source}: IBD spectrum, observed vs fitted", fontsize=14)
                save(fig, folder, f"spectrum_fit_{scenario}_{source}",
                     "Top: one subsample (rep 0), exact Poisson 68% bars, fitted expected counts at the best event order of each graph. "
                     "Bottom: median ± IQR of Pearson residuals across all subsamples.")

            # ---- 4. Ne trajectory ----------------------------------------------
            fig, axes = plt.subplots(2, 4, figsize=(18, 8))
            for ax, branch in zip(axes.flat, BACKBONE):
                ref = reference(c, scenario, branch, collapsed=True)
                t1 = ref["t1"] if ref["t1"] is not None else ref["t0"]+150
                grid = np.linspace(ref["t0"], t1, 400)
                ax.plot(grid, effective_ne(c, scenario, branch, grid, collapsed=True), color="k", lw=2, label="true effective Ne(t)")
                ax.hlines(ref["harmonic"], ref["t0"], t1, colors="k", linestyles=":", lw=1.2, label="harmonic mean")
                ax.hlines(ref["arithmetic"], ref["t0"], t1, colors="k", linestyles="-.", lw=.9, label="arithmetic mean")
                for q, (variant, comp) in enumerate(NE_COMPONENTS):
                    vals = []
                    for r in data[(scenario, source, "omitted_backbone", variant)]:
                        bm = branch_map(r["candidate"])
                        vals += [r[comp]["mean"][j] for j, node in enumerate(r["candidate"]["nodes"]) if bm.get(node) == branch]
                    if vals:
                        qs = np.quantile(vals, [.25, .5, .75])
                        ax.hlines(qs[1], ref["t0"], t1, colors=NE_COLORS[q], lw=2.2, label=NE_LABELS[q].replace("\n", " ")+f" (n={len(vals)})")
                        ax.fill_between([ref["t0"], t1], qs[0], qs[2], color=NE_COLORS[q], alpha=.15)
                ax.set_yscale("log"); ax.set_title(branch); ax.set_xlabel("generations ago"); ax.grid(alpha=.15)
            axes[0, 0].set_ylabel("haploid Ne"); axes[1, 0].set_ylabel("haploid Ne")
            axes[0, 0].legend(fontsize=7, loc="best")
            fig.suptitle(f"{scenario} / {source}: true Ne(t) vs fitted constants on the omitted-loop backbone (g15)", fontsize=14)
            save(fig, folder, f"ne_trajectory_{scenario}_{source}",
                 "Fitted lines = median of posterior means across subsamples, band = IQR, drawn over the TRUE branch interval. "
                 "Collapsed b = leaf growth, then loop 1/Σp²/N, then anc_b.")

        # ---- 5. IBD vs SNP size split ---------------------------------------------
        fig, axes = plt.subplots(2, 4, figsize=(17, 7.5))
        for ax, branch in zip(axes.flat, BACKBONE):
            groups = []
            for source in SOURCES:
                vals = []
                for r in data[(scenario, source, "omitted_backbone", "poisson_separate")]:
                    bm = branch_map(r["candidate"])
                    vals += [np.log(r["Ne_ibd"]["mean"][j]/r["Ne_snp"]["mean"][j]) for j, node in enumerate(r["candidate"]["nodes"]) if bm.get(node) == branch]
                groups.append(vals)
            ref = reference(c, scenario, branch, collapsed=True)
            box(ax, groups, 0.0, ["true IBD", "hap-IBD"], ["#236FA6", "#D87921"])
            ax.set_title(f"{branch}   young {ref['young']:,.0f} / old {ref['old']:,.0f}", fontsize=10)
            ax.set_ylabel("log(Ne_IBD / Ne_SNP)")
        fig.suptitle(f"{scenario}: separate-Ne model, size split by data type (backbone g15)", fontsize=14)
        save(fig, folder, f"ne_ibd_vs_snp_{scenario}",
             "Growth prediction: > 0 on growing branches (IBD sees the recent large end, SNP drift the harmonic mean). Constant control should sit at 0.")

        # ---- 6. residual tilt -----------------------------------------------------
        fig, axes = plt.subplots(2, 6, figsize=(21, 7))
        for i_s, source in enumerate(SOURCES):
            for p, (i, j) in enumerate(PAIRS):
                groups = [[r["residual_tilt"][f"{i},{j}"]["slope"] for r in data[(scenario, source, lab, v)]] for lab, v in COMPARISONS]
                box(axes[i_s, p], groups, 0.0)
                axes[i_s, p].set_title(f"{PAIR_NAMES[p]} / {source}", fontsize=10)
                axes[i_s, p].set_ylabel("Pearson residual slope per cM")
        fig.suptitle(f"{scenario}: residual tilt across segment length", fontsize=14)
        save(fig, folder, f"residual_tilt_{scenario}",
             "Slope of Pearson residual vs bin midpoint. A non-zero slope is the constant-Ne fingerprint of a branch whose size changed.")

        # ---- 7. component likelihood ---------------------------------------------
        n_ibd = len(edges(c))-1; n_ibd *= 6
        fig, axes = plt.subplots(2, 4, figsize=(17, 7.5))
        for i_s, source in enumerate(SOURCES):
            for p, (key, div, title) in enumerate((("lp_ibd", 1, "IBD log-lik"), ("lp_snp", 1, "SNP log-lik"),
                                                   ("chi2_ibd", n_ibd, "IBD chi²/n"), ("chi2_snp", 6, "SNP chi²/n"))):
                groups = [[r["likelihood_summaries"][key]["mean"]/div for r in data[(scenario, source, lab, v)] if key in r.get("likelihood_summaries", {})]
                          for lab, v in COMPARISONS]
                box(axes[i_s, p], groups, 1.0 if key.startswith("chi2") else None)
                axes[i_s, p].set_title(f"{title} / {source}", fontsize=10)
        fig.suptitle(f"{scenario}: composite-likelihood components", fontsize=14)
        save(fig, folder, f"component_likelihood_{scenario}",
             "Levels are comparable within a column only. chi²/n = 1 means the data sit inside their own counting noise.")

        # ---- 8. backbone parameters -----------------------------------------------
        params = ["b_split", "left", "right", "root", "b_fraction"]
        fig, axes = plt.subplots(2, 5, figsize=(18, 8))
        for i_s, source in enumerate(SOURCES):
            for p, param in enumerate(params):
                groups = [[r["semantic_parameters"][param]["mean"] for r in data[(scenario, source, lab, v)] if param in r.get("semantic_parameters", {})]
                          for lab, v in COMPARISONS]
                truth = c["fractions"]["b"] if param == "b_fraction" else c["times"][param]
                box(axes[i_s, p], groups, truth)
                axes[i_s, p].set_title(param.replace("_", " ")+f" (truth {truth:g})")
                axes[i_s, p].set_ylabel(source+(" / fraction" if param == "b_fraction" else " / generations"))
        fig.suptitle(f"{scenario}: event parameters on the correct backbone (descriptive)", fontsize=14)
        save(fig, folder, f"parameters_backbone_{scenario}",
             "Under growth the truth lies outside every candidate's model class; these are not recovery targets.")
