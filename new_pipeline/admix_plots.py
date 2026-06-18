"""
Shared comparison plots for the admix-through-b + ancient-admix experiment
series (true-IBD / Hap-IBD / FastSMC versions).

Two figures are produced by the callers (a third, the relative-error box plot,
lives in relative_error.py):

  1. plot_delbo_boxplot      -- dELBO (T_true - T_alt) box plots with means.
  2. plot_param_comparison   -- parameter recovery, IBD-only vs Mixed columns.

Both read the same per-replicate accumulators the run-scripts already build:
    elbo_results[model_label][cm] = list of tuples (one per replicate), each
        tuple holding the per-topology best ELBO (index 0 == true topology).
    admix_params[model_label][topo_index][cm] = list of dicts (one per rep),
        each mapping a parameter key (e.g. 'a_time', 'a_frac') to its
        posterior-mean point estimate.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def plot_delbo_boxplot(elbo_results, cm_values, model_labels, model_colors,
                       comparisons, outpath, suptitle=""):
    """dELBO box plots (one panel per model x comparison), means marked.

    comparisons : list of (label, idx_a, idx_b); each panel plots
                  ELBO[idx_a] - ELBO[idx_b] across cm as box plots.
    """
    nrows, ncols = len(comparisons), len(model_labels)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(9 * ncols, 6 * nrows), squeeze=False)
    n_cm = len(cm_values)

    for r, (clabel, idx_a, idx_b) in enumerate(comparisons):
        for c, mlabel in enumerate(model_labels):
            ax = axes[r][c]
            color = model_colors[c]

            data, positions = [], []
            for ci, cm in enumerate(cm_values):
                elbos = elbo_results[mlabel][cm]
                if len(elbos) > 0:
                    arr = np.array(elbos)
                    d = arr[:, idx_a] - arr[:, idx_b]
                    d = d[np.isfinite(d)]
                    if d.size:
                        data.append(d)
                        positions.append(ci)

            if data:
                ax.boxplot(
                    data, positions=positions, widths=0.6,
                    showmeans=True, showfliers=False, patch_artist=True,
                    boxprops=dict(facecolor=color, alpha=0.5, edgecolor="k"),
                    medianprops=dict(color="k"),
                    meanprops=dict(marker="D", markerfacecolor="k",
                                   markeredgecolor="k", markersize=5),
                )

            ax.axhline(0, color="k", ls="--", lw=1, alpha=0.5)
            ax.set_xlim(-0.5, n_cm - 0.5)
            ax.set_xticks(range(n_cm))
            ax.set_xticklabels([str(cm) for cm in cm_values], fontsize=10)
            ax.set_xlabel("Genome length (cM)", fontsize=11)
            ax.set_ylabel("delta-ELBO", fontsize=11)
            ax.set_title(f"{mlabel}: {clabel}", fontsize=11, fontweight="bold")
            ax.grid(axis="y", alpha=0.15)

            # % of replicates with dELBO > 0 (true topology preferred)
            ylo = ax.get_ylim()[0]
            for ci, cm in enumerate(cm_values):
                elbos = elbo_results[mlabel][cm]
                if len(elbos) > 0:
                    arr = np.array(elbos)
                    d = arr[:, idx_a] - arr[:, idx_b]
                    d = d[np.isfinite(d)]
                    if d.size:
                        ax.text(ci, ylo, f"{100.0 * np.mean(d > 0):.0f}%",
                                ha="center", va="bottom", fontsize=8,
                                color="green", fontweight="bold")

    if suptitle:
        fig.suptitle(suptitle, fontsize=12, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.93])
    else:
        fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")
    return outpath


def plot_param_comparison(admix_params, cm_values, param_rows, models,
                          topo_colors, outpath, suptitle=""):
    """Parameter recovery: one row per parameter, one column per model.

    param_rows : list of (row_label, param_key, truth, topo_entries) where
                 topo_entries is a list of (topo_index, display_name).
    models     : the model labels to show as columns (e.g. ['IBD-only','Mixed']).
    """
    nrows, ncols = len(param_rows), len(models)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(8 * ncols, 5 * nrows), squeeze=False)

    for r, (rlabel, pkey, truth, topo_entries) in enumerate(param_rows):
        for c, mlabel in enumerate(models):
            ax = axes[r][c]
            n_g = len(topo_entries)
            bar_width = 0.7 / n_g

            for gi, (ti, tdisp) in enumerate(topo_entries):
                color = topo_colors[ti]
                data, positions = [], []
                for ci, cm in enumerate(cm_values):
                    vals = [p[pkey] for p in admix_params[mlabel][ti][cm]
                            if pkey in p]
                    if len(vals) > 0:
                        data.append(np.array(vals))
                        positions.append(ci + (gi - (n_g - 1) / 2) * bar_width)

                if data:
                    vp = ax.violinplot(data, positions=positions,
                                       widths=bar_width * 0.85,
                                       showmeans=True, showextrema=False)
                    for body in vp["bodies"]:
                        body.set_facecolor(color)
                        body.set_alpha(0.55)
                        body.set_edgecolor("k")
                        body.set_linewidth(0.5)
                    vp["cmeans"].set_color(color)
                    vp["cmeans"].set_linewidth(1.5)

            ax.axhline(truth, color="red", ls="--", lw=1.2)
            ax.set_xlim(-0.5, len(cm_values) - 0.5)
            ax.set_xticks(range(len(cm_values)))
            ax.set_xticklabels([str(cm) for cm in cm_values], fontsize=9)
            ax.set_xlabel("Genome length (cM)", fontsize=10)
            ax.set_ylabel(pkey, fontsize=10)
            ax.set_title(f"{mlabel}: {rlabel}", fontsize=10, fontweight="bold")
            ax.grid(axis="y", alpha=0.15)

            if c == 0:
                handles = [Patch(facecolor=topo_colors[ti], alpha=0.55,
                                 label=name) for ti, name in topo_entries]
                handles.append(plt.Line2D([0], [0], color="red", ls="--",
                                          lw=1.2, label="truth"))
                ax.legend(handles=handles, fontsize=8, loc="best")

    if suptitle:
        fig.suptitle(suptitle, fontsize=13, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
    else:
        fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")
    return outpath
