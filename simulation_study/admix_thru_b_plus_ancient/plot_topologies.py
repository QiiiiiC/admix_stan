"""The three candidate graphs, with every event and branch named as in the
parameter figures (`parameters_*.png`, `parameter_estimates.csv`).

Layout rule, the same for every graph: leaves sit at fixed x, a merge's parent
sits at the midpoint of its children, and an admixture's two sources sit at
fixed offsets either side of the admixed branch.  Leaf order is chosen per
graph so no edges cross.  Vertical spacing is schematic (one level per true
event time); events with no counterpart in the truth sit between levels.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from study import CANDIDATES, FRACTIONS, SEMANTICS, build, load_config

INK, MUTED, BRANCH = "#0b0b0b", "#52514e", "#3d3d3a"
ADMIX = "#eb6834"      # admixture source branches (reference categorical slot 2)
# Every label sits on a white backing above the edges, so text is never crossed by a line.
BACK = dict(boxstyle="round,pad=0.15", fc="white", ec="none")
LEVELS = ["loop_split", "loop_donor_merge", "loop_close", "ancient_admix", "left", "right", "root"]

# Per graph: leaf x positions, admixture source offsets, and y level of every event
# (by the node it creates/splits).  Levels 1..7 are the true event times.
LAYOUT = {
    "T_true": dict(leaves={"a": 0.6, "b": 1.6, "c": 3.1, "d": 4.1},
                   admix={"a": (-0.35, 0.45), "c": (-0.45, 0.4)},
                   level={"a": 1, "bP": 2, "ab": 3, "c": 4, "left": 5, "right": 6, "root": 7}),
    "T_alt1": dict(leaves={"a": 0.6, "b": 1.6, "c": 3.1, "d": 4.1},
                   admix={"c": (-0.45, 0.4)},
                   level={"ab": 3, "c": 4, "left": 5, "right": 6, "root": 7}),
    "T_alt2": dict(leaves={"a": 0.6, "b": 1.6, "d": 3.1, "c": 4.1},
                   admix={"a": (-0.35, 0.45)},
                   level={"a": 1, "bP": 2, "ab": 3, "abd": 5.5, "root": 7}),
}
TITLES = {"T_true": "T_true (generating graph)", "T_alt1": "T_alt1: no recent loop",
          "T_alt2": "T_alt2: no ancient admixture"}


def geometry(name):
    """x of every node, and (y_start, y_end) of every branch."""
    d = build(name); lay = LAYOUT[name]
    x = dict(lay["leaves"]); y0 = {leaf: 0.0 for leaf in x}; y1 = {}
    for ev in d.ordered_events:
        if ev["type"] == "MERGE":
            level = lay["level"][ev["parent"]]
            for child in ev["children"]:
                y1[child] = level
            x[ev["parent"]] = sum(x[c] for c in ev["children"]) / 2
            y0[ev["parent"]] = level
        else:
            level = lay["level"][ev["child"]]
            y1[ev["child"]] = level
            for source, dx in zip(ev["parents"], lay["admix"][ev["child"]]):
                x[source] = x[ev["child"]] + dx
                y0[source] = level
    root = next(n for n in d.nodes if n not in y1)
    y1[root] = 7.45
    return d, x, y0, y1


def draw(ax, name, c):
    d, x, y0, y1 = geometry(name)
    names = SEMANTICS[name]
    sources = {p: ev["child"] for ev in d.ordered_events if ev["type"] == "ADMIXTURE" for p in ev["parents"]}
    # branches: vertical run, then the connector to the parent (merge) or to the sources (admixture)
    for node in d.nodes:
        color = ADMIX if node in sources else BRANCH
        ax.plot([x[node], x[node]], [y0[node], y1[node]], color=color, lw=2.2, solid_capstyle="round")
        mid = (y0[node] + min(y1[node], 7.2)) / 2
        ax.text(x[node] + 0.06, mid, node, fontsize=7.5, color=MUTED, va="center", ha="left",
                style="italic", bbox=BACK, zorder=6)
    for ev in d.ordered_events:
        if ev["type"] == "MERGE":
            p = ev["parent"]
            for child in ev["children"]:
                ax.plot([x[child], x[p]], [y0[p], y0[p]], color=BRANCH, lw=2.2, solid_capstyle="round")
            node, level = p, y0[p]
        else:
            child = ev["child"]
            for s in ev["parents"]:
                ax.plot([x[child], x[s]], [y1[child]] * 2, color=ADMIX, lw=2.2, ls=(0, (3, 1.5)))
            node, level = child, y1[child]
            label, key = FRACTIONS[child]
            p1, p2 = ev["parents"]
            has = child in names
            ax.text(x[p1] - 0.05, level + 0.28, label if has else "f", fontsize=7.5, color=ADMIX, ha="right",
                    bbox=BACK, zorder=6)
            ax.text(x[p2] + 0.05, level + 0.28, f"1 − {label}" if has else "1 − f", fontsize=7.5, color=ADMIX,
                    ha="left", bbox=BACK, zorder=6)
        ax.scatter([x[node]], [level], s=28, color=INK, zorder=5)
        if node in names:
            text = f"{names[node]}  ({c['times'][names[node]]:g})" if name == "T_true" else names[node]
            weight = "bold"
        else:
            text, weight = f"{node}: no true counterpart", "normal"
        ha, dx = ("right", -0.1) if ev["type"] == "ADMIXTURE" and x[node] > 2.5 else ("left", 0.1)
        ax.text(x[node] + dx, level - 0.16, text, fontsize=8.5, color=INK, weight=weight, ha=ha, va="top",
                bbox=BACK, zorder=7)
    for leaf, lx in LAYOUT[name]["leaves"].items():
        ax.text(lx, -0.3, leaf, fontsize=14, weight="bold", ha="center", va="top", color=INK)
    ax.set_title(TITLES[name], fontsize=11, color=INK, loc="left")
    ax.set_xlim(-0.2, 4.9); ax.set_ylim(-0.8, 7.7); ax.set_xticks([])
    ax.spines[["top", "right", "bottom"]].set_visible(False)


def topology_figure(c, destination):
    destination = Path(destination); destination.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(17, 7.5), sharey=True)
    for ax, name in zip(axes, CANDIDATES):
        draw(ax, name, c)
    axes[0].set_yticks(range(8), ["0"] + [f"{c['times'][k]:g}" for k in LEVELS])
    axes[0].set_ylabel("generations before present (schematic spacing)")
    f = c["fractions"]
    fig.text(0.01, 0.015,
             f"Bold: event times as named in parameters_*.png (T_true shows the truth in brackets).  "
             f"Orange dashed: admixture; fractions go to the FIRST source (loop_fraction = {f['loop']:g}, "
             f"ancient_fraction = {f['ancient']:g}).  Grey italic: branch (node) names.  "
             f"Ne: one shared size on every branch (truth {c['haploid_ne']:,.0f} haploid).",
             fontsize=8.5, color=MUTED)
    fig.suptitle("Candidate graphs", x=0.01, ha="left", fontsize=13, color=INK)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    for ext in ("png", "pdf"):
        fig.savefig(destination / f"topology_candidates.{ext}", dpi=200 if ext == "png" else None,
                    bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    topology_figure(load_config(), here / "summary" / "figures")
