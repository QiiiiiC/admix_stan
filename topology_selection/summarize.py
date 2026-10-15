"""report.md and figures for runs/<out>.

Per dataset: the graph posterior from the full enumeration (softmax of best-path ELBO
plus log prior), and for every search how many fits it needed to find the best graph
and how much of the posterior mass the graphs it scored cover.

  search.png     share of searches that have found the best graph, by number of fits
  landscape.png  the ground truth: every graph's ELBO relative to the best, per dataset
  graphs.png     the graphs that matter (top five per dataset, the walks' traps, T_null)
"""
from __future__ import annotations

from collections import Counter
import itertools
import json

import numpy as np
from scipy.special import logsumexp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import graphs as gr
from methods.demography import DemographicTopology

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
STYLE = {  # method -> (colour, line style, marker): identity never by colour alone
    "walk, beta=1": ("#2a6fdb", "-", "o"),
    "walk, beta=0.3": ("#1a9e77", "--", "s"),
    "walk, beta=0.1": ("#7b5ea7", ":", "^"),
    "walk + relocate, beta=1": ("#c2185b", "-", "P"),
    "walk + relocate, beta=0.3": ("#c2185b", "--", "X"),
    "greedy": ("#e08a1e", "-", "D"),
    "greedy + relocate": ("#e08a1e", "--", "d"),
    "stepwise (trees, then all)": ("#6b6b6b", "-.", "v"),
}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "axes.axisbelow": True})


def named_graphs():
    t = DemographicTopology(list("abcd"))
    t.add_admixture_event("b", "b1", "b2"); t.add_merge_event("a", "b1", "ab"); t.add_merge_event("b2", "c", "cb")
    t.add_merge_event("cb", "d", "cbd"); t.add_merge_event("ab", "cbd", "root")
    n = DemographicTopology(list("abcd"))
    n.add_merge_event("a", "b", "ab"); n.add_merge_event("c", "d", "cd"); n.add_merge_event("ab", "cd", "root")
    return {gr.from_topology(t).key(): "T_true", gr.from_topology(n).key(): "T_null"}


def segments(g, start, node):
    """Each non-root branch as a vertical (x, y0, y1) and a horizontal (y, x0, x1) piece."""
    con, root, out = g.consumer(), g.root(), []
    for b, (x0, y0) in start.items():
        if b != root:
            x1, y1 = node[con[b]]
            out.append((b, (x0, y0, y1), (y1, min(x0, x1), max(x0, x1))))
    return out


def clutter(g, start, node):
    """Edge crossings, plus near-coincident vertical lines (they read as one line)."""
    seg, bad = segments(g, start, node), 0
    for b, v, _ in seg:
        for b2, _, h in seg:
            if b != b2 and h[1] + 1e-9 < v[0] < h[2] - 1e-9 and v[1] + 1e-9 < h[0] < v[2] - 1e-9:
                bad += 1
    for i, (_, v, _) in enumerate(seg):
        for _, v2, _ in seg[i + 1:]:
            if abs(v[0] - v2[0]) < 0.3 and min(v[2], v2[2]) - max(v[1], v2[1]) > 0.2:
                bad += 1
    return bad


def layout(g):
    """Schematic coordinates: an event's level is 1 + the highest level below it.  Leaf
    order and admixture sides are chosen to minimise crossings and overlapping lines,
    then horizontal travel.  Returns the start point of every branch and the point of
    every event."""
    cre, best = g.creator(), None
    level = {}
    for e in g.listing():
        level[e] = 1 + max([level[cre[b]] for b in g.events[e][1] if b in cre], default=0)
    admixes = [e for e, (k, _, _) in g.events.items() if k == "admix"]
    for perm in itertools.permutations(g.leaves):
        for flips in itertools.product((1, -1), repeat=len(admixes)):
            side = dict(zip(admixes, flips))
            start = {b: (float(perm.index(b)), 0.0) for b in g.leaves}
            node, travel = {}, 0.0
            for e in g.listing():
                kind, ins, outs = g.events[e]
                y = float(level[e])
                if kind == "merge":
                    x = float(np.mean([start[b][0] for b in ins]))
                    travel += sum(abs(start[b][0] - x) for b in ins)
                    start[outs[0]] = (x, y)
                else:
                    x = start[ins[0]][0]
                    start[outs[0]] = (x - 0.3 * side[e], y); start[outs[1]] = (x + 0.3 * side[e], y)
                node[e] = (x, y)
            cost = (clutter(g, start, node), travel)
            if best is None or cost < best[0]:
                best = (cost, start, node)
    return best[1], best[2]


def draw_graph(g, ax):
    """Branches as elbows (up from their start, across to the event that consumes them);
    admixture sources dashed until they merge; admixture events as diamonds."""
    start, node = layout(g)
    con, root = g.consumer(), g.root()
    source = {b for k, _, outs in g.events.values() if k == "admix" for b in outs}
    for b, (x0, y0) in start.items():
        if b == root:
            ax.plot([x0, x0], [y0, y0 + 0.6], color=INK, lw=2, solid_capstyle="round")
            continue
        x1, y1 = node[con[b]]
        ls = (0, (3, 2)) if b in source else "-"
        colour = "#c2185b" if b in source else INK
        ax.plot([x0, x0, x1], [y0, y1, y1], color=colour, lw=2, ls=ls, solid_capstyle="round")
    for e, (x, y) in node.items():
        kind, ins, outs = g.events[e]
        if kind == "admix":
            ax.plot([start[outs[0]][0], start[outs[1]][0]], [y, y], color="#c2185b", lw=2)
            ax.plot(x, y, marker="D", ms=7, color="#c2185b", mec="white", mew=1.5, zorder=3)
        else:
            ax.plot(x, y, marker="o", ms=4, color=INK, zorder=3)
    for b in g.leaves:
        ax.text(start[b][0], -0.35, b, ha="center", va="top", fontsize=11, color=INK)
    top = max(y for _, y in node.values())
    ax.set_xlim(-0.6, len(g.leaves) - 0.4); ax.set_ylim(-1.0, top + 0.8)
    ax.axis("off")


def short(ds):
    return f"{ds[0][0]}{ds[1]}"


def landscape(out, datasets, name, seen, scores, gallery_ids):
    fig, axes = plt.subplots(1, len(datasets), figsize=(4.2 * len(datasets), 3.8), squeeze=False, sharey=True)
    for ax, ds in zip(axes[0], datasets):
        score = scores[ds]
        order = sorted(seen, key=lambda k: -score[k])
        top = score[order[0]]
        for k_admix, (colour, marker, lab) in {1: ("#2a6fdb", "o", "one admixture"),
                                               0: ("#e08a1e", "^", "tree")}.items():
            idx = [i for i, k in enumerate(order) if seen[k].n_admix() == k_admix]
            ax.scatter(np.array(idx) + 1, [score[order[i]] - top for i in idx], s=18, color=colour, marker=marker,
                       edgecolor="white", linewidth=0.5, label=lab, zorder=2)
        for i, k in enumerate(order):
            if k in gallery_ids and (i < 3 or gallery_ids[k] in ("T_true", "T_null")):
                ax.annotate(gallery_ids[k], (i + 1, score[k] - top), xytext=(5, -3 if i else 4),
                            textcoords="offset points", fontsize=8, color=INK)
        ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1)
        ax.set_xlabel("rank of the graph (of 153)"); ax.set_title(name(ds), loc="left", fontsize=10)
        ax.set_ylim(-2e4, 3)
    axes[0, 0].set_ylabel("ELBO minus the best graph's (nats)")
    handles, labs = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labs, frameon=False, fontsize=9, loc="lower center", ncol=2)
    fig.suptitle("Ground truth: every graph in the space, fitted once", x=0.01, ha="left", fontsize=12)
    fig.tight_layout(rect=(0, 0.07, 1, 0.94))
    fig.savefig(out / "landscape.png", dpi=170, bbox_inches="tight"); plt.close(fig)


def gallery(out, datasets, seen, scores, ids):
    keys = list(ids)
    ncol = 4
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.3 * nrow), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, k in zip(axes.ravel(), keys):
        draw_graph(seen[k], ax)
        gaps = "  ".join(f"{short(ds)} {scores[ds][k] - max(scores[ds].values()):+.0f}" for ds in datasets)
        ax.set_title(f"{ids[k]}\n{gaps}", loc="left", fontsize=8.5, color=INK)
    fig.suptitle("The graphs that matter.  Below each name: ELBO minus the best graph's, per dataset "
                 "(g/v = growth/varying design, 100/1000 cM).  Dashed: an admixture's two source lineages.",
                 x=0.01, ha="left", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out / "graphs.png", dpi=170, bbox_inches="tight"); plt.close(fig)


def gallery_ids(datasets, seen, scores, nbrs_score_optima, labels):
    """Top five per dataset, every local optimum of the plain moves, T_null; numbered by mean rank."""
    chosen = set(labels)
    for ds in datasets:
        chosen |= set(sorted(seen, key=lambda k: -scores[ds][k])[:5])
        chosen |= set(nbrs_score_optima[ds])
    rank = {k: np.mean([sorted(seen, key=lambda q: -scores[ds][q]).index(k) for ds in datasets]) for k in chosen}
    ids, n = {}, 0
    for k in sorted(chosen, key=lambda q: (rank[q], q)):     # key breaks ties: same numbering in every run
        if k in labels:
            ids[k] = labels[k]
        else:
            n += 1; ids[k] = f"G{n}"
    return ids


def curves(runs, best, post, grid):
    """Per run: has the best graph been scored by fit n; posterior mass of the scored set."""
    found = np.array([[best in r[:n] for n in grid] for r in runs], float)
    mass = np.array([[sum(post[k] for k in set(r[:n])) for n in grid] for r in runs])
    return found.mean(0), np.median(mass, 0)


def build(out, datasets, name, seen, prior, load_scores, nbrs=None):
    labels = named_graphs()
    scores = {ds: load_scores(out, ds)[0] for ds in datasets}
    optima = {ds: [k for k in seen if all(scores[ds][q] <= scores[ds][k] for q in nbrs[k])] for ds in datasets}
    ids = gallery_ids(datasets, seen, scores, optima, labels)
    landscape(out, datasets, name, seen, scores, ids)
    gallery(out, datasets, seen, scores, ids)
    lines = ["# Topology selection on the four-leaf data", "",
             f"Space: every graph on a, b, c, d with at most one admixture, no bubble, no ghost source: "
             f"{len(seen)} graphs ({sum(g.n_admix() == 0 for g in seen.values())} trees).  Each fitted once with the "
             "per-branch model, Pathfinder 3 starts x 1 path; score = best-path ELBO.  Prior: uniform over the number of "
             "admixtures, then uniform within it.  Posterior = softmax(score + log prior) over the whole space.", "",
             "Searches are replayed on the cached scores; **fits** = distinct graphs a search scored (each would be one "
             "Pathfinder fit in a live run).  Walks start from a random tree; greedy and stepwise start from each of the "
             "15 trees.", ""]
    fig, axes = plt.subplots(2, len(datasets), figsize=(4.2 * len(datasets), 7), squeeze=False)
    grid = np.arange(1, len(seen) + 1)
    for col, ds in enumerate(datasets):
        score, rows = load_scores(out, ds)
        keys = list(seen)
        lp = np.array([score[k] + prior[k] for k in keys])
        post = dict(zip(keys, np.exp(lp - logsumexp(lp))))
        order = sorted(keys, key=lambda k: -score[k])
        best = max(keys, key=lambda k: post[k])
        secs = [r["seconds"] for r in rows]
        failed = sum(any("error" in s for s in r["starts"]) for r in rows)
        lines += [f"## {name(ds)}", "",
                  f"{len(rows)} graphs fitted: median {np.median(secs):.1f} s per graph (wall, 10 in parallel), "
                  f"{failed} with a failed start.", "",
                  "| rank | graph | admixtures | ELBO | Δ to best | posterior |", "|---:|---|---:|---:|---:|---:|"]
        for i, k in enumerate(order[:8]):
            tag = f" **{labels[k]}**" if k in labels else ""
            lines.append(f"| {i + 1} | `{k}`{tag} | {seen[k].n_admix()} | {score[k]:.1f} | "
                         f"{score[k] - score[order[0]]:+.1f} | {post[k]:.3g} |")
        for k, lab in labels.items():
            if k not in order[:8]:
                lines.append(f"| {order.index(k) + 1} | `{k}` **{lab}** | {seen[k].n_admix()} | {score[k]:.1f} | "
                             f"{score[k] - score[order[0]]:+.1f} | {post[k]:.3g} |")
        runs = json.loads((out / name(ds) / "search.json").read_text())
        lines += ["", "A search selects the highest-scoring graph it scored (names as in graphs.png).", "",
                  "| search | runs | found the best graph | fits to find it (median, 90%) | fits used (median) | "
                  "posterior mass scored at the end (median) | selected graph (share of runs) |",
                  "|---|---:|---:|---:|---:|---:|---|"]
        for method, rs in runs.items():
            hit = [r.index(best) + 1 for r in rs if best in r]
            q = f"{np.median(hit):.0f}, {np.percentile(hit, 90):.0f}" if hit else "–"
            mass = [sum(post[k] for k in set(r)) for r in rs]
            picked = Counter(ids.get(max(r, key=score.get), "other") for r in rs)
            lines.append(f"| {method} | {len(rs)} | {len(hit) / len(rs):.0%} | {q} | {np.median([len(r) for r in rs]):.0f} | "
                         f"{np.median(mass):.3f} | " + ", ".join(f"{g} {n / len(rs):.0%}" for g, n in picked.most_common()) + " |")
            found, mass_c = curves(rs, best, post, grid)
            colour, ls, mk = STYLE[method]
            kw = dict(color=colour, ls=ls, lw=2, marker=mk, ms=4, markevery=max(1, len(grid) // 8), label=method)
            axes[0, col].plot(grid, found, **kw)
            axes[1, col].plot(grid, mass_c, **kw)
        lines.append("")
        axes[0, col].set_title(name(ds), loc="left", fontsize=10)
        axes[0, col].set_ylabel("runs that found the best graph")
        axes[1, col].set_ylabel("posterior mass scored (median)")
        for ax in axes[:, col]:
            ax.set_xscale("log"); ax.set_xlabel(f"fits (distinct graphs scored; space = {len(seen)})")
            ax.set_ylim(-0.03, 1.03)
    handles, names = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, names, frameon=False, fontsize=9, loc="lower center", ncol=4)
    fig.suptitle("Graph search on cached Pathfinder ELBOs: four-leaf space, at most one admixture",
                 x=0.01, ha="left", fontsize=12)
    fig.tight_layout(rect=(0, 0.07, 1, 0.96))
    fig.savefig(out / "search.png", dpi=170, bbox_inches="tight"); plt.close(fig)
    (out / "report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
