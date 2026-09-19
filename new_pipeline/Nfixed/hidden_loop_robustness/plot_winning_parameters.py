"""Boxplots of fitted posterior means across subsamples, for selected Poisson graphs."""
import argparse
import csv
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".build" / "matplotlib"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COLORS = ["#2676B8", "#E58A25"]
SCENARIOS = ["hidden_loop", "no_loop"]
SOURCES = ["true", "hapibd"]
VARIANTS = ["poisson_shared", "poisson_separate"]
NODES = ["a", "b", "c", "b.1", "b.2", "n1", "n2", "root"]
NODE_LABELS = ["a", "b", "c", "b.1", "b.2", "Ancestor of a + b.1", "Ancestor of c + b.2", "Root"]


def collect(run):
    with (run / "summary" / "topology_rankings.csv").open() as f:
        winners = [r for r in csv.DictReader(f) if r["rank"] == "1" and r["variant"] in VARIANTS]
    records = []
    for row in winners:
        folder = (run / row["scenario"] / f"pool_{int(row['pool']):03d}" / "fits"
                  / f"rep_{int(row['replicate']):03d}" / row["source"] / row["variant"])
        orders = [json.loads(p.read_text()) for p in folder.glob("*/result.json")]
        orders = [r for r in orders if r["status"] == "ok" and r["candidate"]["graph"] == int(row["graph"])]
        fit = max(orders, key=lambda r: r["logz"])
        # This run selects the same graph and event order in all 80 comparisons.
        # Fail explicitly if a future run needs a different parameter mapping.
        if fit["candidate"]["newick"] != "((a,b.1),(b.2,c))" or not fit["truth"]["correct_order"]:
            raise ValueError("Selected graph/order differs: map comparable events explicitly before plotting")
        identity = {k: row[k] for k in ("scenario", "pool", "replicate", "source", "variant")}
        identity.update(candidate=fit["candidate"]["id"], reliable=fit["reliable"])
        values = dict(zip(["b_time", "left_time", "right_time", "root_time"], fit["cumulative_times"]["mean"]))
        values["fraction"] = fit["admixture_fractions"]["mean"][0]
        for name, value in values.items():
            records.append(dict(identity, parameter=name, component="demography", estimate=value))
        for variable, component in (("Ne", "shared"), ("Ne_ibd", "IBD"), ("Ne_snp", "SNP")):
            if variable in fit:
                for node, value in zip(fit["candidate"]["nodes"], fit[variable]["mean"]):
                    records.append(dict(identity, parameter=node, component=component, estimate=value))
    for scenario in SCENARIOS:
        for source in SOURCES:
            for variant in VARIANTS:
                subset = [r for r in records if r["scenario"] == scenario and r["source"] == source
                          and r["variant"] == variant and r["parameter"] == "b_time"]
                if len(subset) != 10 or len({(r["pool"], r["replicate"]) for r in subset}) != 10:
                    raise ValueError(f"Expected 10 unique subsamples for {scenario}/{source}/{variant}")
    return records


def draw_boxes(ax, groups, colors, labels, truth):
    positions = np.arange(1, len(groups)+1)
    bp = ax.boxplot(groups, positions=positions, widths=0.58, patch_artist=True,
                    showfliers=False, medianprops=dict(color="#17212B", linewidth=1.6),
                    whiskerprops=dict(color="#555555"), capprops=dict(color="#555555"))
    for box, color in zip(bp["boxes"], colors):
        box.set_facecolor(color); box.set_alpha(0.38); box.set_edgecolor(color)
    rng = np.random.default_rng(20260909)
    for pos, values, color in zip(positions, groups, colors):
        ax.scatter(pos+rng.uniform(-0.17, 0.17, len(values)), values, s=23,
                   color=color, edgecolors="white", linewidths=0.35, zorder=3)
    ax.axhline(truth, color="#A52A35", linestyle="--", linewidth=1.5, zorder=2)
    ax.set_xticks(positions, labels, fontsize=8)
    ax.grid(axis="y", alpha=0.17); ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)


def save(fig, out, stem):
    fig.savefig(out / (stem+".png"), dpi=190)
    fig.savefig(out / (stem+".pdf"))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=HERE / "runs" / "default")
    args = parser.parse_args(); run = args.run.resolve()
    records = collect(run)
    config = json.loads((run / "manifest.json").read_text())["config"]
    out = run / "summary" / "parameter_boxplots"; out.mkdir(exist_ok=True)
    with (out / "plotted_estimates.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, list(records[0])); writer.writeheader(); writer.writerows(records)
    def values(scenario, source, variant, parameter, component):
        return [r["estimate"] for r in records if r["scenario"] == scenario and r["source"] == source
                and r["variant"] == variant and r["parameter"] == parameter and r["component"] == component]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titleweight": "semibold"})
    params = ["b_time", "left_time", "right_time", "root_time", "fraction"]
    titles = ["b admixture time", "a + b.1 merge time", "c + b.2 merge time", "Root merge time", "b ancestry fraction toward a"]
    truths = [config["times"][k] for k in ("b_split", "left", "right", "root")] + [config["fractions"]["b"]]
    fig, axes = plt.subplots(2, 5, figsize=(19, 8.3))
    for i, scenario in enumerate(SCENARIOS):
        for j, (param, title, truth) in enumerate(zip(params, titles, truths)):
            groups = [values(scenario, src, variant, param, "demography") for src in SOURCES for variant in VARIANTS]
            draw_boxes(axes[i,j], groups, [COLORS[0]]*2+[COLORS[1]]*2,
                       ["Shared\nTrue IBD", "Separate\nTrue IBD", "Shared\nhap-IBD", "Separate\nhap-IBD"], truth)
            axes[i,j].set_title(title+f"\nTruth = {truth:g}", fontsize=11)
            axes[i,j].set_ylabel(("Hidden loop\n" if i == 0 else "No-loop control\n")+
                                ("Generations ago" if j < 4 else "Ancestry fraction"))
            if param == "fraction":
                axes[i,j].set_ylim(0, 1)
    # Identical y-scales between scenarios for each parameter.
    for j in range(4):
        limits = [ax.get_ylim() for ax in axes[:,j]]
        for ax in axes[:,j]: ax.set_ylim(min(x[0] for x in limits), max(x[1] for x in limits))
    fig.suptitle("Winning topology: ((a, b.1), (b.2, c)) | Poisson shared vs separate Ne", fontsize=17, y=0.98)
    fig.text(0.5, 0.025, "Each point = one subsample's fitted posterior mean (10 per box). Dashed red = simulated truth.\n"
             "Boxes show median and interquartile range; whiskers = 1.5 IQR. Subsamples overlap; inference diagnostics failed, so estimates are provisional.",
             ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.09, 1, 0.93)); save(fig, out, "demographic_parameters")
    for scenario in SCENARIOS:
        fig, axes = plt.subplots(2, 4, figsize=(18, 9))
        for ax, node, label in zip(axes.flat, NODES, NODE_LABELS):
            groups = [values(scenario, src, variant, node, component) for src in SOURCES
                      for variant, component in (("poisson_shared", "shared"), ("poisson_separate", "IBD"), ("poisson_separate", "SNP"))]
            draw_boxes(ax, groups, [COLORS[0]]*3+[COLORS[1]]*3,
                       ["Shared", "Sep.\nIBD", "Sep.\nSNP"]*2, config["haploid_ne"])
            ax.set_title(label); ax.set_ylabel("Haploid Ne")
            ax.ticklabel_format(axis="y", style="plain", useOffset=False)
            ax.set_ylim(bottom=0)
        fig.suptitle("Branch Ne in the winning topology — "+("hidden c loop" if scenario == "hidden_loop" else "no-loop control"), fontsize=17)
        fig.legend(handles=[Line2D([0],[0], marker="s", color=COLORS[0], linestyle="", label="True IBD input (left three boxes)"),
                            Line2D([0],[0], marker="s", color=COLORS[1], linestyle="", label="hap-IBD input (right three boxes)"),
                            Line2D([0],[0], color="#A52A35", linestyle="--", label="Truth: Ne = 15,000")],
                   loc="upper center", bbox_to_anchor=(0.5,0.945), ncol=3, frameon=False)
        fig.text(0.5,0.025,"Each point = one subsample's fitted posterior mean (10 per box). Sep. IBD / SNP = separate-Ne model components.\n"
                 "Boxes are across subsamples, not posterior intervals. Subsamples overlap; numerical reliability criteria were not met.",ha="center",fontsize=10)
        fig.tight_layout(rect=(0,0.085,1,0.90)); save(fig,out,"branch_ne_"+scenario)
    print(f"Saved 3 PNGs, 3 PDFs and {len(records)} plotted estimates to {out}")


if __name__ == "__main__":
    main()
