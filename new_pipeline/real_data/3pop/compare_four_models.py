"""Compare simple- and grid-Ne variants across three-population topologies."""
import argparse
import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr


VARIANTS = [
    "poisson_shared_ne",
    "poisson_separate_ne",
    "normal_shared_ne",
    "normal_separate_ne",
    "poisson_grid_shared_ne",
    "poisson_grid_separate_ne",
    "normal_grid_shared_ne",
    "normal_grid_separate_ne",
]
SIMPLE_VARIANTS = VARIANTS[:4]

LABELS = {
    "poisson_shared_ne": "Poisson, shared Ne",
    "poisson_separate_ne": "Poisson, separate Ne",
    "normal_shared_ne": "Normal, shared Ne",
    "normal_separate_ne": "Normal, separate Ne",
    "poisson_grid_shared_ne": "Poisson, shared Ne, grid",
    "poisson_grid_separate_ne": "Poisson, separate Ne, grid",
    "normal_grid_shared_ne": "Normal, shared Ne, grid",
    "normal_grid_separate_ne": "Normal, separate Ne, grid",
}


def softmax(x):
    z = np.asarray(x, float) - np.max(x)
    w = np.exp(z)
    return w / w.sum()


def admixture_fraction(row):
    fractions = row.get("admixture_fractions", [])
    return float(fractions[0]) if fractions else None


def format_fraction(row, digits=3):
    fraction = admixture_fraction(row)
    return "-" if fraction is None else f"{fraction:.{digits}f}"


def closest_non_admix_tree(row, pops, all_rows, threshold=0.01):
    """Collapse a boundary admixture edge and identify its displayed tree."""
    fraction = admixture_fraction(row)
    if fraction is None or min(fraction, 1.0 - fraction) >= threshold:
        return None

    tree_by_sisters = {}
    for candidate in all_rows:
        if candidate["n_admix"] != 0:
            continue
        first_merge = next(e for e in candidate["events"] if e["type"] == "MERGE")
        tree_by_sisters[frozenset(first_merge["children"])] = candidate

    admixture = next(e for e in row["events"] if e["type"] == "ADMIXTURE")
    retained = admixture["parents"][0 if fraction > 0.5 else 1]
    states = {pop: {pop} for pop in pops}
    states[admixture["child"]] = set()
    for parent in admixture["parents"]:
        states[parent] = {admixture["child"]} if parent == retained else set()

    for event in row["events"]:
        if event["type"] != "MERGE":
            continue
        left = states.get(event["children"][0], set())
        right = states.get(event["children"][1], set())
        merged = left | right
        states[event["parent"]] = merged
        if left and right and len(merged) == 2:
            return tree_by_sisters.get(frozenset(merged))
    return None


def format_closest_tree(row, pops, all_rows):
    tree = closest_non_admix_tree(row, pops, all_rows)
    return "-" if tree is None else f"[{tree['index']:02d}] {tree['newick']}"


def load_variant(root, key):
    path = os.path.join(root, "comparison", key, "all_fits.json")
    with open(path) as fh:
        blob = json.load(fh)
    rows = sorted(blob["results"], key=lambda r: -r["elbo"])
    top = rows[0]["elbo"]
    probs = softmax([r["elbo"] for r in rows])
    for rank, (row, prob) in enumerate(zip(rows, probs), 1):
        row["rank"] = rank
        row["d_elbo"] = row["elbo"] - top
        row["p_model"] = float(prob)
    return blob["pops"], rows


def write_variant_outputs(root, key, pops, rows):
    out = os.path.join(root, "comparison", key)
    os.makedirs(out, exist_ok=True)
    txt = [f"TOPOLOGY RANKING: {LABELS[key]} | {'/'.join(pops)}", ""]
    txt.append(f"{'rank':>4} {'idx':>3} {'topology':<35} {'f':>8} "
               f"{'closest tree when boundary':<25} {'ELBO':>12} "
               f"{'dELBO':>10} {'ELBO-w':>8} {'X2 IBD':>9} {'X2 SNP':>9}")
    for r in rows:
        ci = r["chi2_ibd"] / max(r.get("n_ibd_obs", 1), 1)
        cs = r["chi2_snp"] / max(r.get("n_snp_obs", 1), 1)
        txt.append(f"{r['rank']:>4} {r['index']:>3} {r['newick']:<35} "
                   f"{format_fraction(r, 5):>8} "
                   f"{format_closest_tree(r, pops, rows):<25} "
                   f"{r['elbo']:>12.2f} {r['d_elbo']:>10.2f} "
                   f"{r['p_model']:>8.4f} {ci:>9.2f} {cs:>9.2f}")
    with open(os.path.join(out, "ranking.txt"), "w") as fh:
        fh.write("\n".join(txt) + "\n")

    with open(os.path.join(out, "elbo_table.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "index", "newick", "n_admix", "admixture_fraction",
                    "closest_tree_index", "closest_tree_newick",
                    "elbo", "d_elbo",
                    "elbo_softmax_weight", "elbo_se", "ess", "chi2_ibd_per_n",
                    "chi2_snp_per_n", "mode", "map_start", "map_lp",
                    "n_map_success", "n_pathfinder_modes", "elbo_mode_range",
                    "seed", "secs"])
        for r in rows:
            closest = closest_non_admix_tree(r, pops, rows)
            w.writerow([
                r["rank"], r["index"], r["newick"], r["n_admix"],
                "" if admixture_fraction(r) is None else f"{admixture_fraction(r):.8f}",
                "" if closest is None else closest["index"],
                "" if closest is None else closest["newick"],
                f"{r['elbo']:.6f}", f"{r['d_elbo']:.6f}",
                f"{r['p_model']:.8g}", f"{r['elbo_se']:.6f}",
                f"{r['ess']:.3f}",
                f"{r['chi2_ibd']/max(r.get('n_ibd_obs',1),1):.6f}",
                f"{r['chi2_snp']/max(r.get('n_snp_obs',1),1):.6f}",
                r.get("mode", ""), r.get("map_start", ""),
                "" if "map_lp" not in r else f"{r['map_lp']:.6f}",
                r.get("n_map_success", ""), r.get("n_pathfinder_modes", ""),
                "" if "elbo_mode_min" not in r else
                f"{r['elbo_mode_max'] - r['elbo_mode_min']:.6f}",
                r["seed"], f"{r['secs']:.1f}",
            ])

    y = np.arange(len(rows))[::-1]
    colors = ["#2f6f9f" if r["n_admix"] == 0 else "#b5483f" for r in rows]
    fig, ax = plt.subplots(figsize=(17, 7.5))
    ax.barh(y, [r["d_elbo"] for r in rows], color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels([f"[{r['index']:02d}] {r['newick']}" for r in rows], fontsize=8)
    for yi, r in zip(y, rows):
        ax.text(1.015, yi, format_fraction(r, 5), transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=8, family="monospace")
        ax.text(1.11, yi, format_closest_tree(r, pops, rows),
                transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=8)
    ax.text(1.015, 1.012, "f", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9, fontweight="bold")
    ax.text(1.11, 1.012, "closest non-admix tree", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=9, fontweight="bold")
    ax.axvline(0, color="black", lw=.8)
    ax.set_xlabel("ELBO - best within model (nats; natural-log units, 0 = best)")
    ax.set_title(f"{LABELS[key]} | {' / '.join(pops)}")
    fig.tight_layout(rect=(0, 0, 0.72, 1))
    fig.savefig(os.path.join(out, "elbo_ranking.png"), dpi=140)
    plt.close(fig)

    report = [f"# {LABELS[key]} topology ranking", "",
              f"Populations: **{' / '.join(pops)}**", "",
              "The weight is a softmax of Pathfinder ELBOs, not a posterior model probability.",
              "Each topology is screened from dispersed MAP starts; distinct high-MAP modes "
              "are then fitted independently with Pathfinder. Search diagnostics are retained "
              "in `fit.json` and `elbo_table.csv`.",
              "`f` is the fitted ancestry fraction from the first named source branch; "
              "tree topologies have no fraction. The closest tree is shown when "
              "`min(f, 1-f) < 0.01`.", "",
              "| rank | topology | f | closest non-admix tree | ELBO | dELBO | ELBO weight |",
              "|---|---|---:|---|---:|---:|---:|"]
    for r in rows:
        link = f"../../{r['name']}/{key}/report.md"
        closest = format_closest_tree(r, pops, rows)
        report.append(f"| {r['rank']} | [`{r['newick']}`]({link}) | {format_fraction(r, 5)} | "
                      f"`{closest}` | "
                      f"{r['elbo']:+.1f} | {r['d_elbo']:+.1f} | {r['p_model']:.4f} |")
    report += ["", "![ranking](elbo_ranking.png)", ""]
    with open(os.path.join(out, "report.md"), "w") as fh:
        fh.write("\n".join(report))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--simple-only", action="store_true",
                    help="Regenerate only the four completed non-grid variants.")
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.join(here, args.tag)

    by_variant = {}
    pops = None
    variants = SIMPLE_VARIANTS if args.simple_only else VARIANTS
    for key in variants:
        p, rows = load_variant(root, key)
        if pops is None:
            pops = p
        elif p != pops:
            raise ValueError(f"population order differs for {key}: {p} != {pops}")
        by_variant[key] = rows
        write_variant_outputs(root, key, pops, rows)

    by_index = {key: {r["index"]: r for r in rows}
                for key, rows in by_variant.items()}
    indices = sorted(set.intersection(*(set(v) for v in by_index.values())))
    comparison = os.path.join(root, "comparison")

    if args.simple_only:
        with open(os.path.join(comparison, "four_model_topology_table.csv"),
                  "w", newline="") as fh:
            w = csv.writer(fh)
            header = ["index", "newick"]
            for key in variants:
                header += [f"rank_{key}", f"f_{key}",
                           f"closest_tree_{key}", f"elbo_{key}",
                           f"d_elbo_{key}"]
            w.writerow(header)
            for idx in indices:
                base = by_index[variants[0]][idx]
                row = [idx, base["newick"]]
                for key in variants:
                    r = by_index[key][idx]
                    fraction = admixture_fraction(r)
                    closest = closest_non_admix_tree(r, pops, by_variant[key])
                    row += [r["rank"],
                            "" if fraction is None else f"{fraction:.8f}",
                            "" if closest is None else closest["newick"],
                            f"{r['elbo']:.6f}", f"{r['d_elbo']:.6f}"]
                w.writerow(row)

        x = np.arange(len(indices))
        fig, axes = plt.subplots(1, 2, figsize=(15, 5), sharex=True, sharey=True)
        colors = {"poisson_shared_ne": "#2166ac",
                  "poisson_separate_ne": "#b2182b",
                  "normal_shared_ne": "#1b7837",
                  "normal_separate_ne": "#762a83"}
        for ax, likelihood in zip(axes, ("poisson", "normal")):
            for key in [k for k in variants if k.startswith(likelihood)]:
                ax.plot(x, [by_index[key][i]["rank"] for i in indices], "o-",
                        ms=4, lw=1.2, color=colors[key], label=LABELS[key])
            ax.invert_yaxis()
            ax.set_xticks(x)
            ax.set_xticklabels([f"{i:02d}" for i in indices])
            ax.set_xlabel("topology index (categorical, not ELBO)")
            ax.set_title(f"{likelihood.title()} topology ranks")
            ax.legend(fontsize=8)
        axes[0].set_ylabel("rank within model (1 = best)")
        fig.suptitle(f"Four completed simple-Ne models | {' / '.join(pops)}")
        fig.tight_layout()
        fig.savefig(os.path.join(comparison, "four_model_comparison.png"), dpi=140)
        plt.close(fig)

        report = [f"# Four completed simple-Ne models on {' / '.join(pops)}", "",
                  "These are the completed multistart reruns. Grid-Ne results are "
                  "excluded until that batch finishes.", "",
                  "## Best topology within each model", "",
                  "| model | winner | ELBO (nats) | runner-up dELBO | mode range | ESS |",
                  "|---|---|---:|---:|---:|---:|"]
        for key in variants:
            rows = by_variant[key]
            winner = rows[0]
            mode_range = (winner.get("elbo_mode_max", winner["elbo"])
                          - winner.get("elbo_mode_min", winner["elbo"]))
            report.append(f"| {LABELS[key]} | `{winner['newick']}` | "
                          f"{winner['elbo']:+.1f} | {rows[1]['d_elbo']:+.1f} | "
                          f"{mode_range:.1f} | {winner['ess']:.1f} |")
        report += ["", "## Rank consistency", "",
                   "The horizontal position in the figure is only the topology index. "
                   "ELBO ranking is on the vertical axis; individual ranking plots use "
                   "horizontal `ELBO - best` in natural-log units (nats).", "",
                   "| comparison | Spearman rho | top-5 overlap |",
                   "|---|---:|---:|"]
        pairs = [("Poisson shared vs separate", variants[0], variants[1]),
                 ("Normal shared vs separate", variants[2], variants[3]),
                 ("Shared Ne: Poisson vs Normal", variants[0], variants[2]),
                 ("Separate Ne: Poisson vs Normal", variants[1], variants[3])]
        for label, left, right in pairs:
            lr = {r["index"]: r["rank"] for r in by_variant[left]}
            rr = {r["index"]: r["rank"] for r in by_variant[right]}
            rho = spearmanr([lr[i] for i in indices],
                            [rr[i] for i in indices]).statistic
            overlap = len(set(sorted(lr, key=lr.get)[:5]) &
                          set(sorted(rr, key=rr.get)[:5]))
            report.append(f"| {label} | {rho:+.3f} | {overlap}/5 |")
        report += ["", "![four-model comparison](four_model_comparison.png)", ""]
        with open(os.path.join(comparison, "four_model_report.md"), "w") as fh:
            fh.write("\n".join(report))
        return

    with open(os.path.join(comparison, "eight_model_topology_table.csv"),
              "w", newline="") as fh:
        w = csv.writer(fh)
        header = ["index", "newick"]
        for key in VARIANTS:
            header += [f"rank_{key}", f"f_{key}", f"closest_tree_{key}",
                       f"elbo_{key}", f"d_elbo_{key}"]
        w.writerow(header)
        for idx in indices:
            base = by_index[VARIANTS[0]][idx]
            row = [idx, base["newick"]]
            for key in VARIANTS:
                r = by_index[key][idx]
                fraction = admixture_fraction(r)
                closest = closest_non_admix_tree(r, pops, by_variant[key])
                row += [r["rank"], "" if fraction is None else f"{fraction:.8f}",
                        "" if closest is None else closest["newick"],
                        f"{r['elbo']:.6f}", f"{r['d_elbo']:.6f}"]
            w.writerow(row)

    fig, ax = plt.subplots(2, 2, figsize=(15, 10), sharex=True)
    x = np.arange(len(indices))
    styles = {
        "poisson_shared_ne": ("#2166ac", "o", "-"),
        "poisson_separate_ne": ("#b2182b", "s", "-"),
        "poisson_grid_shared_ne": ("#67a9cf", "^", "--"),
        "poisson_grid_separate_ne": ("#ef8a62", "D", "--"),
        "normal_shared_ne": ("#1b7837", "o", "-"),
        "normal_separate_ne": ("#762a83", "s", "-"),
        "normal_grid_shared_ne": ("#7fbf7b", "^", "--"),
        "normal_grid_separate_ne": ("#af8dc3", "D", "--"),
    }
    likelihood_keys = {
        "Poisson": [k for k in VARIANTS if k.startswith("poisson")],
        "Normal": [k for k in VARIANTS if k.startswith("normal")],
    }
    for col, (likelihood, keys) in enumerate(likelihood_keys.items()):
        for key in keys:
            color, marker, line = styles[key]
            ax[0, col].plot(
                x, [by_index[key][i]["rank"] for i in indices], line,
                marker=marker, ms=4, lw=1.2, color=color, label=LABELS[key])
        ax[0, col].set_yticks([1, 5, 10, 15, 21])
        ax[0, col].invert_yaxis()
        ax[0, col].set_title(f"{likelihood} topology ranks")
        ax[0, col].legend(fontsize=8)
    ax[0, 0].set_ylabel("rank within variant (1 = best)")

    gain_pairs = [
        ("Poisson shared", "poisson_shared_ne", "poisson_grid_shared_ne", "#2166ac"),
        ("Poisson separate", "poisson_separate_ne", "poisson_grid_separate_ne", "#b2182b"),
        ("Normal shared", "normal_shared_ne", "normal_grid_shared_ne", "#1b7837"),
        ("Normal separate", "normal_separate_ne", "normal_grid_separate_ne", "#762a83"),
    ]
    for col, likelihood in enumerate(("Poisson", "Normal")):
        for label, simple, grid, color in gain_pairs:
            if not label.startswith(likelihood):
                continue
            gain = [by_index[grid][i]["elbo"] - by_index[simple][i]["elbo"]
                    for i in indices]
            ax[1, col].plot(x, gain, "o-", ms=4, lw=1.2, color=color, label=label)
        ax[1, col].axhline(0, color="black", lw=.8)
        ax[1, col].set_title(f"{likelihood}: grid minus simple Ne")
        ax[1, col].set_xticks(x)
        ax[1, col].set_xticklabels([f"{i:02d}" for i in indices])
        ax[1, col].set_xlabel("topology index")
        ax[1, col].legend(fontsize=8)
    ax[1, 0].set_ylabel("ELBO gain (nats)")
    fig.suptitle(f"Eight-model comparison | {' / '.join(pops)}")
    fig.tight_layout()
    fig.savefig(os.path.join(comparison, "eight_model_comparison.png"), dpi=140)
    plt.close(fig)

    report = [f"# Eight-model comparison on {' / '.join(pops)}", "",
              "## Best topology within each variant", "",
              "| variant | topology | ELBO | dELBO runner-up |",
              "|---|---|---:|---:|"]
    for key in VARIANTS:
        rows = by_variant[key]
        report.append(f"| {LABELS[key]} | `{rows[0]['newick']}` | "
                      f"{rows[0]['elbo']:+.1f} | {rows[1]['d_elbo']:+.1f} |")
    report += ["", "## Ne-model evidence gains", "",
               "These differences are valid only within the same likelihood.", "",
               "| likelihood | comparison | best-to-best gain | same winner? |",
               "|---|---|---:|---:|"]
    evidence_pairs = [
        ("Poisson", "simple: separate - shared", "poisson_shared_ne", "poisson_separate_ne"),
        ("Poisson", "grid: separate - shared", "poisson_grid_shared_ne", "poisson_grid_separate_ne"),
        ("Poisson", "shared: grid - simple", "poisson_shared_ne", "poisson_grid_shared_ne"),
        ("Poisson", "separate: grid - simple", "poisson_separate_ne", "poisson_grid_separate_ne"),
        ("Normal", "simple: separate - shared", "normal_shared_ne", "normal_separate_ne"),
        ("Normal", "grid: separate - shared", "normal_grid_shared_ne", "normal_grid_separate_ne"),
        ("Normal", "shared: grid - simple", "normal_shared_ne", "normal_grid_shared_ne"),
        ("Normal", "separate: grid - simple", "normal_separate_ne", "normal_grid_separate_ne"),
    ]
    for likelihood, label, base_key, expanded_key in evidence_pairs:
        base = by_variant[base_key][0]
        expanded = by_variant[expanded_key][0]
        gain = expanded["elbo"] - base["elbo"]
        same = "yes" if expanded["index"] == base["index"] else "no"
        report.append(f"| {likelihood} | {label} | {gain:+.1f} | {same} |")
    report += ["", "![comparison](eight_model_comparison.png)", "",
               "## Rank consistency", "",
               "ELBO values are comparable between shared- and separate-Ne models "
               "within one likelihood, but not between Poisson and Normal because "
               "they are masses/densities for different observed summaries.", "",
               "| comparison | Spearman rank correlation | top-5 overlap |",
               "|---|---:|---:|"]
    rank_pairs = [
        ("Poisson simple: shared vs separate", "poisson_shared_ne", "poisson_separate_ne"),
        ("Poisson grid: shared vs separate", "poisson_grid_shared_ne", "poisson_grid_separate_ne"),
        ("Poisson shared: simple vs grid", "poisson_shared_ne", "poisson_grid_shared_ne"),
        ("Poisson separate: simple vs grid", "poisson_separate_ne", "poisson_grid_separate_ne"),
        ("Normal simple: shared vs separate", "normal_shared_ne", "normal_separate_ne"),
        ("Normal grid: shared vs separate", "normal_grid_shared_ne", "normal_grid_separate_ne"),
        ("Normal shared: simple vs grid", "normal_shared_ne", "normal_grid_shared_ne"),
        ("Normal separate: simple vs grid", "normal_separate_ne", "normal_grid_separate_ne"),
        ("Simple shared: Poisson vs Normal", "poisson_shared_ne", "normal_shared_ne"),
        ("Simple separate: Poisson vs Normal", "poisson_separate_ne", "normal_separate_ne"),
        ("Grid shared: Poisson vs Normal", "poisson_grid_shared_ne", "normal_grid_shared_ne"),
        ("Grid separate: Poisson vs Normal", "poisson_grid_separate_ne", "normal_grid_separate_ne"),
    ]
    for label, left, right in rank_pairs:
        lr = {r["index"]: r["rank"] for r in by_variant[left]}
        rr = {r["index"]: r["rank"] for r in by_variant[right]}
        rho = spearmanr([lr[i] for i in indices], [rr[i] for i in indices]).statistic
        overlap = len(set(sorted(lr, key=lr.get)[:5]) & set(sorted(rr, key=rr.get)[:5]))
        report.append(f"| {label} | {rho:+.3f} | {overlap}/5 |")

    report += ["", "## Winner diagnostics", "",
               "| variant | IBD chi2/n | SNP chi2/n | admixture fraction | interpretation |",
               "|---|---:|---:|---:|---|"]
    for key in VARIANTS:
        winner = by_variant[key][0]
        f = winner.get("admixture_fractions", [])
        if not f:
            fraction, interpretation = "-", "tree"
        else:
            fraction = f"{f[0]:.4f}"
            interpretation = ("collapsed/near-tree" if min(f[0], 1.0 - f[0]) < 0.05
                              else "non-boundary admixture")
        report.append(
            f"| {LABELS[key]} | "
            f"{winner['chi2_ibd']/max(winner.get('n_ibd_obs', 1), 1):.2f} | "
            f"{winner['chi2_snp']/max(winner.get('n_snp_obs', 1), 1):.2f} | "
            f"{fraction} | {interpretation} |")

    recent_rows = []
    for key in [k for k in VARIANTS if "_grid_" in k]:
        winner = by_variant[key][0]
        components = ([('shared', 'Ne_recent')] if not winner.get('separate_ne') else
                      [('IBD', 'Ne_recent_ibd'), ('SNP', 'Ne_recent_snp')])
        for component, field in components:
            recent = np.asarray(winner[field], float)
            for pop, values in zip(pops, recent):
                recent_rows.append({
                    "variant": key,
                    "topology_index": winner["index"],
                    "component": component,
                    "population": pop,
                    "ne_0_5": values[0],
                    "ne_5_10": values[1],
                    "recent_ratio": values[0] / values[1],
                    "first_event_generation": winner["times"][0],
                })
    with open(os.path.join(comparison, "recent_ne_winners.csv"),
              "w", newline="") as fh:
        fields = list(recent_rows[0])
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(recent_rows)

    report += ["", "## Recent Ne in grid-model winners", "",
               "The ratio is Ne(0-5 generations) / Ne(5-10 generations); values above "
               "one indicate growth toward the present.", "",
               "| variant | component | population | Ne 0-5 | Ne 5-10 | ratio | first event |",
               "|---|---|---|---:|---:|---:|---:|"]
    for row in recent_rows:
        report.append(
            f"| {LABELS[row['variant']]} | {row['component']} | {row['population']} | "
            f"{row['ne_0_5']:.0f} | {row['ne_5_10']:.0f} | "
            f"{row['recent_ratio']:.2f} | {row['first_event_generation']:.1f} |")

    count = np.asarray(by_variant["poisson_grid_shared_ne"][0]["ibd_count"])
    n_zero = sum(np.count_nonzero(count[:, i, j] == 0)
                 for i in range(count.shape[1]) for j in range(i, count.shape[2]))
    n_unique = count.shape[0] * count.shape[1] * (count.shape[1] + 1) // 2
    normal_winner = by_variant["normal_grid_separate_ne"][0]
    normal_floor = np.asarray(normal_winner.get("ibd_theory_se_pred", []))
    n_floor = (sum(np.count_nonzero(normal_floor[:, i, j] <= 1.0001e-12)
                   for i in range(normal_floor.shape[1])
                   for j in range(i, normal_floor.shape[2]))
               if normal_floor.size else 0)
    poisson_winner = by_variant["poisson_grid_separate_ne"][0]
    pf = poisson_winner.get("admixture_fractions", [np.nan])[0]
    winner_lines = []
    collapsed_indices = []
    for key in VARIANTS:
        winner = by_variant[key][0]
        closest = closest_non_admix_tree(winner, pops, by_variant[key])
        collapsed_indices.append(None if closest is None else closest["index"])
        winner_lines.append(
            f"{LABELS[key]} selects topology {winner['index']} "
            f"(`{winner['newick']}`), f={format_fraction(winner, 4)}"
            + (f", collapsing to tree {closest['index']} (`{closest['newick']}`)."
               if closest is not None else "."))
    common_collapsed = (collapsed_indices[0] if collapsed_indices and
                        all(i == collapsed_indices[0] for i in collapsed_indices)
                        else None)
    report += ["", "## Interpretation", "",
               "Poisson and Normal ELBO levels are not subtracted from each other: the two "
               "likelihoods are densities/masses for different summaries and therefore use "
               "different base measures and units. Compare their topology ranks, residuals, "
               "and shared-versus-separate-Ne gains instead.", "",
               f"The data contain **{n_zero}/{n_unique} empty unique pair-by-bin cells**. "
               f"The winning grid separate-Ne Normal fit places {n_floor}/{n_unique} unique cells at its "
               "hard `1e-12` theory-SE floor. Consequently, its all-bin CLT likelihood is "
               "being used most aggressively exactly where the CLT is least justified.", "",
               f"The Poisson winner has admixture fraction **{pf:.4f}**. It is therefore "
               "best read as a near-tree model with an extra branch breakpoint if the fraction "
               "remains near a boundary, not as evidence for substantial admixture.", "",
               "The grid comparison directly tests whether 0-10 generation growth explains the "
               "topology preference. Rank agreement and the winner-fraction diagnostics above "
               "should be considered together; a high ELBO for boundary admixture is evidence "
               "for remaining Ne misspecification rather than for gene flow.", "",
               "## Conclusion", "",
               " ".join(winner_lines), "",
               (f"All eight winners collapse to the same non-admixture topology "
                f"**{common_collapsed}**. The raw graph labels differ between shared- and "
                "separate-Ne parameterizations, but their boundary fractions make them "
                "equivalent at the population-tree level."
                if common_collapsed is not None else
                "The winning graphs do not all collapse to one non-admixture topology; "
                "their fractions and collapsed-tree labels should therefore be compared "
                "individually."), "",
               f"The grid separate-Ne Poisson winner has IBD chi2/n "
               f"**{poisson_winner['chi2_ibd']/max(poisson_winner.get('n_ibd_obs', 1), 1):.2f}** "
               f"but fraction **{pf:.4f}**. Its good count calibration does not turn that "
               "boundary edge into admixture evidence; the graph is acting as a tree with an "
               "additional branch-specific Ne change.", "",
               f"The grid separate-Ne Normal winner has fraction "
               f"**{normal_winner.get('admixture_fractions', [float('nan')])[0]:.4f}**, but its "
               f"IBD chi2/n is **{normal_winner['chi2_ibd']/max(normal_winner.get('n_ibd_obs', 1), 1):.2f}** "
               f"and {n_floor}/{n_unique} cells use the variance floor. Its topology result is "
               "therefore sensitivity evidence, not a reliable resolution of the Poisson result.", "",
               "The next misspecification test should use fixed absolute Ne breakpoints beyond "
               "generation 10 and truncate each branch trajectory at its event time. That lets "
               "events occur inside the grid, avoiding an artificial lower bound on the first "
               "event while testing whether the near-tree admixture edge disappears.", ""]
    with open(os.path.join(comparison, "report.md"), "w") as fh:
        fh.write("\n".join(report))


if __name__ == "__main__":
    main()
