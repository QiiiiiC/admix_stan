"""Draw the configured generating histories; vertical spacing is schematic."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def topology_figures(config, destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    t = config["times"]
    for loop in (False, True):
        fig, ax = plt.subplots(figsize=(10, 8))
        blue, orange = "#236FA6", "#D87921"

        def edge(points, color=blue):
            ax.plot(*zip(*points), color=color, lw=2.5, solid_capstyle="round")

        def node(x, y, label, dx=0.12, dy=0.08):
            ax.scatter([x], [y], s=35, color=blue, zorder=5)
            ax.text(x+dx, y+dy, label, fontsize=11, va="bottom")

        edge([(0, 0), (0, 4), (.65, 4)])
        edge([(4, 0), (4, 5), (3.35, 5)])
        edge([(2, 3), (1.25, 3), (1.25, 4), (.65, 4)])
        edge([(2, 3), (2.75, 3), (2.75, 5), (3.35, 5)])
        edge([(.65, 4), (.65, 6), (2, 6)])
        edge([(3.35, 5), (3.35, 6), (2, 6)])
        edge([(2, 6), (2, 6.45)])
        if loop:
            edge([(2, 0), (2, 1)])
            edge([(2, 1), (1.45, 1), (1.45, 2), (2, 2)], orange)
            edge([(2, 1), (2.55, 1), (2.55, 2), (2, 2)], orange)
            edge([(2, 2), (2, 3)])
            node(2, 1, "Loop opens", dx=.65, dy=-.12)
            node(2, 2, "Loop closes: anc_b", dx=.25)
            f = config["fractions"]["loop"]
            ax.text(1.25, 1.5, f"loop1\n{f:.0%}", ha="right", va="center", color=orange)
            ax.text(2.75, 1.5, f"loop2\n{1-f:.0%}", ha="left", va="center", color=orange)
        else:
            edge([(2, 0), (2, 3)])
        node(2, 3, "b admixture", dx=.05, dy=-.35)
        node(.65, 4, "a + b.1", dx=-.15)
        node(3.35, 5, "b.2 + c", dx=.12)
        node(2, 6, "Root", dx=.15)
        f = config["fractions"]["b"]
        ax.text(1.23, 3.25, f"b.1: {f:.0%}", ha="right", color=blue, fontsize=11)
        ax.text(3.05, 3.85, f"b.2: {1-f:.0%}", color=blue, fontsize=11)
        for x, label in [(0, "a"), (2, "b"), (4, "c")]:
            ax.scatter([x], [0], s=60, color=blue)
            ax.text(x, -.25, label, ha="center", fontsize=16, weight="bold")
        ys = [0, 1, 2, 3, 4, 5, 6] if loop else [0, 3, 4, 5, 6]
        ts = [0, t["loop_open"], t["loop_close"], t["b_split"], t["left"], t["right"], t["root"]]
        ax.set_yticks(ys, [str(ts[y]) for y in ys])
        ax.set_ylabel("Generations before present (schematic spacing)", fontsize=11)
        ax.set_xticks([])
        ax.set_xlim(-.65, 4.85)
        ax.set_ylim(-.5, 6.6)
        ax.grid(axis="y", alpha=.15)
        ax.spines[["top", "right", "bottom"]].set_visible(False)
        title = "b-loop topology — graph 22, order 1" if loop else "Backbone topology — graph 15, order 1"
        ax.set_title(title, fontsize=17, pad=18)
        fig.text(.5, .035,
                 "Configured generating values, not fitted estimates. Time runs upward into the past.\n"
                 f"All branches: haploid Ne = {config['haploid_ne']:,.0f}. Fractions label backward ancestry allocation.",
                 ha="center", fontsize=10)
        fig.tight_layout(rect=(0, .09, 1, 1))
        name = "topology_b_loop" if loop else "topology_backbone"
        for ext in ("png", "pdf"):
            fig.savefig(destination / f"{name}.{ext}", dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    run = here / "runs" / "default"
    topology_figures(json.loads((run / "manifest.json").read_text())["config"], run / "summary" / "figures")
