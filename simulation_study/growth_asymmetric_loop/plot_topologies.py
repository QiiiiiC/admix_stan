"""Draw the generating history for each scenario with branch sizes; spacing is schematic."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from study import branch_sizes, SCENARIOS


def topology_figures(config, destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    t = config["times"]
    for scenario in SCENARIOS:
        sizes = branch_sizes(config, scenario)
        fig, ax = plt.subplots(figsize=(11, 8.5))
        blue, orange, red = "#236FA6", "#D87921", "#A42E32"

        def edge(points, color=blue, lw=2.5):
            ax.plot(*zip(*points), color=color, lw=lw, solid_capstyle="round")

        def node(x, y, label, dx=0.12, dy=0.08):
            ax.scatter([x], [y], s=35, color=blue, zorder=5)
            ax.text(x+dx, y+dy, label, fontsize=11, va="bottom")

        def size(x, y, name, ha="left", color=red):
            young, old = sizes[name]
            txt = f"{young:,.0f}" if young == old else f"{young:,.0f} → {old:,.0f}"
            ax.text(x, y, txt, ha=ha, va="center", fontsize=9, color=color,
                    bbox=dict(boxstyle="round,pad=.15", fc="white", ec="none", alpha=.8))

        edge([(0, 0), (0, 4), (.65, 4)]);  size(.12, 2.0, "a")
        edge([(4, 0), (4, 5), (3.35, 5)]); size(3.88, 2.5, "c", ha="right")
        edge([(2, 3), (1.25, 3), (1.25, 4), (.65, 4)]); size(1.13, 3.5, "b1", ha="right")
        edge([(2, 3), (2.75, 3), (2.75, 5), (3.35, 5)]); size(2.87, 4.0, "b2")
        edge([(.65, 4), (.65, 6), (2, 6)]); size(.77, 5.0, "left")
        edge([(3.35, 5), (3.35, 6), (2, 6)]); size(3.23, 5.5, "right", ha="right")
        edge([(2, 6), (2, 6.45)]); size(2.12, 6.3, "root")
        edge([(2, 0), (2, 1)]); size(2.12, .5, "b")
        edge([(2, 1), (1.45, 1), (1.45, 2), (2, 2)], orange); size(1.33, 1.5, "loop1", ha="right")
        edge([(2, 1), (2.55, 1), (2.55, 2), (2, 2)], orange); size(2.67, 1.5, "loop2")
        edge([(2, 2), (2, 3)]); size(2.12, 2.5, "anc_b")
        node(2, 1, "loop opens", dx=.65, dy=-.12)
        node(2, 2, "loop closes: anc_b", dx=.25)
        f = config["fractions"]["loop"]
        ax.text(1.45, .78, f"{f:.0%}", ha="center", va="top", color=orange, fontsize=9)
        ax.text(2.55, .78, f"{1-f:.0%}", ha="center", va="top", color=orange, fontsize=9)
        node(2, 3, "b admixture", dx=.05, dy=-.35)
        node(.65, 4, "a + b.1", dx=-.15)
        node(3.35, 5, "b.2 + c", dx=.12)
        node(2, 6, "root", dx=.15)
        f = config["fractions"]["b"]
        ax.text(1.23, 2.78, f"b.1: {f:.0%}", ha="right", color=blue, fontsize=10)
        ax.text(2.77, 2.78, f"b.2: {1-f:.0%}", color=blue, fontsize=10)
        for x, label in [(0, "a"), (2, "b"), (4, "c")]:
            ax.scatter([x], [0], s=60, color=blue)
            ax.text(x, -.25, label, ha="center", fontsize=16, weight="bold")
        ys = [0, 1, 2, 3, 4, 5, 6]
        ts = [0, t["loop_open"], t["loop_close"], t["b_split"], t["left"], t["right"], t["root"]]
        ax.set_yticks(ys, [str(v) for v in ts])
        ax.set_ylabel("generations before present (schematic spacing)", fontsize=11)
        ax.set_xticks([]); ax.set_xlim(-.65, 4.85); ax.set_ylim(-.5, 6.6)
        ax.grid(axis="y", alpha=.15); ax.spines[["top", "right", "bottom"]].set_visible(False)
        ax.set_title(f"generating history — {scenario}", fontsize=17, pad=18)
        note = ("Red labels: haploid Ne at the young end → old end of each branch (exponential in between).\n"
                "a grows toward the present across its whole leaf branch; b across its leaf branch back to the loop."
                if scenario == "growth" else
                "Red labels: haploid Ne, constant on every branch. Same topology and times as `growth`.")
        fig.text(.5, .03, note, ha="center", fontsize=10)
        fig.tight_layout(rect=(0, .08, 1, 1))
        for ext in ("png", "pdf"):
            fig.savefig(destination / f"topology_generating_{scenario}.{ext}", dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    run = here / "runs" / "default"
    topology_figures(json.loads((run / "manifest.json").read_text())["config"], run / "summary" / "figures")
