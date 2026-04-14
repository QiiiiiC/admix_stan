"""
Topology swap experiment (Nvarying):
  Can the three models (IBD, SNP, Mixed) distinguish a true pure tree
  from a wrong admixture topology?

Data simulated under T1 (true):
  Pure tree ((a,b),c),d)
  Merge times: a+b=20, ab+c=700, abc+d=1700
  Ne(a)=200 (haploid), all other nodes=20000 (haploid)

T2 (wrong — admixture topology):
  d is treated as admixed: d -> dPa1 (from a lineage) + dPa2 (from b lineage)
  a + dPa1 -> aDmix   (first merge)
  b + dPa2 -> bDmix   (second merge)
  aDmix + bDmix -> subAnc
  subAnc + c -> anc

Plot: delta-ELBO(T1 - T2) and w(T1) for IBD-only, SNP-only, and Mixed.
"""

import sys, os
_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _PARENT)
_MODELS = os.path.join(_PARENT, "models")

import numpy as np
import matplotlib.pyplot as plt
import re

from simulation_methods import (
    simulate_msprime,
    build_ibd_stan_data,
    build_snp_stan_data,
    build_mixed_stan_data,
)
from inference_methods import (
    prepare_snp_blocks,
    pool_snp_blocks,
    resample_snp_covariance,
)
from ibd_jackknife import (
    calculate_ibd_blocks_mrca,
    pool_multiple_simulations,
    resample_ibd_with_jackknife_variance,
)
from demography import DemographicTopology
from cmdstanpy import CmdStanModel


def extract_elbo(fit):
    """Extract the best ELBO from pathfinder stdout files."""
    best_elbo = -np.inf
    for stdout_file in fit._runset.stdout_files:
        try:
            with open(stdout_file, 'r') as f:
                for line in f:
                    m = re.search(r'Best Iter:.*?ELBO \(([-+\d.eE]+)\)', line)
                    if m:
                        elbo_val = float(m.group(1))
                        if elbo_val > best_elbo:
                            best_elbo = elbo_val
        except (OSError, ValueError):
            continue
    return best_elbo if np.isfinite(best_elbo) else np.nan


# ====================================================================
# 0. Configuration
# ====================================================================
N_REPLICATES = 200
BLOCK_SIZE_CM = 50.0
CM_PER_UNIT = 1e-4
RECOMB_RATE = 1e-6
MUT_RATE = 1e-6
SAMPLES_PER_POP = {'a': 15, 'b': 15, 'c': 15, 'd': 15}

N_SIMS = 50
SIM_CM_EACH = 50
SIM_SEQ_LEN = SIM_CM_EACH / CM_PER_UNIT

cm_values = [50, 100, 150, 200, 250, 300, 400, 500]

bins = [
    [0.5, 0.55], [0.55, 0.6], [0.6, 0.65], [0.65, 0.7],
    [0.7, 0.8], [0.8, 0.9], [0.9, 1.0],
    [1.0, 1.5], [1.5, 2.0], [2.0, 5.0], [5.0, 8.0],
    [8.0, 20.0], [20.0, BLOCK_SIZE_CM]
]

NE_A = 200       # haploid Ne for node a
NE_REST = 20000  # haploid Ne for all other nodes

SNP_CUTOFF_TIME = 1700
SNP_MIN_MAF = 0.05
FIXED_NE = NE_REST  # passed to build_snp_stan_data (ignored by Nvarying model)


# ====================================================================
# 1. Define topologies
# ====================================================================
print("=" * 60)
print("Defining topologies")
print("=" * 60)

# --- T1 (true): pure tree ((a,b),c),d) ---
dem_sim = DemographicTopology(['a', 'b', 'c', 'd'])
dem_sim.add_merge_event('a', 'b', 'ab')
dem_sim.add_merge_event('ab', 'c', 'abc')
dem_sim.add_merge_event('abc', 'd', 'root')

dem_sim.set_node_ne('a', NE_A // 2)  # msprime uses diploid Ne
for node in ['b', 'c', 'd', 'ab', 'abc', 'root']:
    dem_sim.set_node_ne(node, NE_REST // 2)

dem_sim.set_merge_time('ab', 20)
dem_sim.set_merge_time('abc', 700)
dem_sim.set_merge_time('root', 1700)
dem_sim.finalize_root()

dem_T1 = dem_sim

# --- T2 (wrong): admixture topology with same leaves {a, b, c, d} ---
dem_T2 = DemographicTopology(['a', 'b', 'c', 'd'])
dem_T2.add_admixture_event('b', 'bPa1', 'bPa2')
dem_T2.add_merge_event('a', 'bPa1', 'aBmix')       
dem_T2.add_merge_event('c', 'bPa2', 'cBmix')       
dem_T2.add_merge_event('aBmix', 'cBmix', 'subAnc')
dem_T2.add_merge_event('subAnc', 'd', 'anc')
dem_T2.finalize_root()

topology_dems = [dem_T1, dem_T2]
topology_labels = ["T1 (true: pure tree)", "T2 (wrong: admixture)"]

for lbl, d in zip(topology_labels, topology_dems):
    print(f"  {lbl}: n_events={len(d.ordered_events)}, n_admix={d.n_admix}, "
          f"n_nodes={len(d.nodes)}")


# ====================================================================
# 2. Simulate and pool blocks (using true topology)
# ====================================================================
packed_list = []
n_blocks_list = []
snp_blocks_list = []
pair_info = None
pop_ids = None

for sim_i in range(N_SIMS):
    print(f"\n--- Simulation {sim_i + 1}/{N_SIMS} ({SIM_CM_EACH} cM) ---")
    mts = simulate_msprime(
        dem_sim,
        sequence_length=SIM_SEQ_LEN,
        recombination_rate=RECOMB_RATE,
        mutation_rate=MUT_RATE,
        samples_per_pop=SAMPLES_PER_POP,
        seed=42 + sim_i,
    )

    packed, n_blocks, pop_samples, pids, pi = calculate_ibd_blocks_mrca(
        mts,
        bins=bins,
        block_size_cm=BLOCK_SIZE_CM,
        cm_per_unit=CM_PER_UNIT,
    )
    packed_list.append(packed)
    n_blocks_list.append(n_blocks)
    if pair_info is None:
        pair_info = pi
        pop_ids = pids

    snp_blks, _ = prepare_snp_blocks(
        mts, dem_sim,
        block_size_cm=BLOCK_SIZE_CM,
        cm_per_unit=CM_PER_UNIT,
        cutoff_time=SNP_CUTOFF_TIME,
        min_maf=SNP_MIN_MAF,
    )
    snp_blocks_list.append(snp_blks)

pooled_ibd, total_blocks = pool_multiple_simulations(packed_list, n_blocks_list)
pooled_snp = pool_snp_blocks(snp_blocks_list)

print(f"\nTotal pool: {total_blocks} blocks of {BLOCK_SIZE_CM} cM "
      f"= {total_blocks * BLOCK_SIZE_CM} cM from {N_SIMS} independent simulations")


# ====================================================================
# 3. Compile Stan models (Nvarying)
# ====================================================================
ibd_stan = CmdStanModel(stan_file=os.path.join(_MODELS, "ibd_model_Nvarying.stan"))
snp_stan = CmdStanModel(stan_file=os.path.join(_MODELS, "snp_model_Nvarying.stan"))
mixed_stan = CmdStanModel(stan_file=os.path.join(_MODELS, "mixed_model_Nvarying.stan"))

model_labels = ["IBD-only", "SNP-only", "Mixed"]
model_stans = [ibd_stan, snp_stan, mixed_stan]
model_colors = ["#378ADD", "#4CAF50", "#E8A838"]


# ====================================================================
# 4. Run replicates: fit T1 vs T2 with all three models
# ====================================================================
# elbo_results[model_label][cm_val] = list of (elbo_T1, elbo_T2) tuples
elbo_results = {label: {cm: [] for cm in cm_values} for label in model_labels}

rng = np.random.default_rng(seed=2025)

for cm_val in cm_values:
    print(f"\n{'='*60}")
    print(f"  cm = {cm_val} cM  ({N_REPLICATES} reps x 2 topologies x 3 models)")
    print(f"{'='*60}")

    n_blocks_needed = max(1, int(round(cm_val / BLOCK_SIZE_CM)))

    for rep in range(N_REPLICATES):
        if (rep + 1) % 10 == 0:
            print(f"  replicate {rep + 1}/{N_REPLICATES}")

        chosen_blocks = rng.choice(total_blocks, size=n_blocks_needed, replace=False)

        # --- Resample data (shared across topologies and models) ---
        ibd_mean, ibd_var = resample_ibd_with_jackknife_variance(
            pooled_ibd,
            n_blocks_total=total_blocks,
            pop_ids=pop_ids,
            pair_info=pair_info,
            bins=bins,
            target_cm=cm_val,
            block_size_cm=BLOCK_SIZE_CM,
            chosen_blocks=chosen_blocks,
        )

        leaves = dem_sim.initial_leaves
        n_haploid = np.array([SAMPLES_PER_POP[p] * 2 for p in leaves])
        w_hat, w_se = resample_snp_covariance(
            pooled_snp, chosen_blocks,
            n_haploid_per_pop=n_haploid,
            se_block_size=50,
        )

        # --- Fit each model x each topology ---
        for m_idx, (m_label, m_stan) in enumerate(zip(model_labels, model_stans)):
            rep_elbos = [np.nan, np.nan]

            for ti, dem_infer in enumerate(topology_dems):
                n_events_t = len(dem_infer.ordered_events)
                n_nodes_t = len(dem_infer.nodes)
                n_admix_t = dem_infer.n_admix

                # Build init dict
                init = {
                    "times": [100.0] * n_events_t,
                    "Ne": [10000.0] * n_nodes_t,
                }
                if n_admix_t > 0:
                    init["admixture_fractions"] = [0.5] * n_admix_t
                if m_label == "Mixed":
                    init["kappa_snp"] = 1.0

                try:
                    if m_label == "IBD-only":
                        sd = build_ibd_stan_data(
                            dem_infer, ibd_mean, ibd_var, bins,
                            T_max=100000, cm=cm_val,
                        )
                    elif m_label == "SNP-only":
                        sd = build_snp_stan_data(
                            dem_infer, w_hat, w_se, effective_N=FIXED_NE,
                        )
                    else:  # Mixed
                        sd = build_mixed_stan_data(
                            dem_infer, ibd_mean, ibd_var, bins, w_hat, w_se,
                            T_max=100000, cm=cm_val,
                        )

                    fit = m_stan.pathfinder(
                        data=sd, inits=init, show_console=False,
                        psis_resample=True,
                    )
                    rep_elbos[ti] = extract_elbo(fit)
                except Exception as e:
                    if rep < 3:
                        print(f"    [WARN] {m_label} {topology_labels[ti]} "
                              f"rep {rep+1}: {e}")

            if all(np.isfinite(rep_elbos)):
                elbo_results[m_label][cm_val].append(tuple(rep_elbos))


# ====================================================================
# 5. Plot: 2x3 grid  (row = metric, col = model)
# ====================================================================
from matplotlib.patches import Patch

fig, axes = plt.subplots(2, 3, figsize=(28, 16))

x = np.arange(len(cm_values))

for m_idx, m_label in enumerate(model_labels):
    color = model_colors[m_idx]
    ax_elbo = axes[0, m_idx]
    ax_w = axes[1, m_idx]

    diff_data = []
    weight_data = []
    for cm_val in cm_values:
        elbos = elbo_results[m_label][cm_val]
        if len(elbos) > 0:
            arr = np.array(elbos)  # (n_reps, 2)
            diff_data.append(arr[:, 0] - arr[:, 1])
            # softmax weight for T1
            max_e = np.max(arr, axis=1, keepdims=True)
            exp_e = np.exp(arr - max_e)
            w = exp_e / exp_e.sum(axis=1, keepdims=True)
            weight_data.append(w[:, 0])
        else:
            diff_data.append(np.array([np.nan]))
            weight_data.append(np.array([np.nan]))

    # --- Top: delta-ELBO ---
    vp1 = ax_elbo.violinplot(
        diff_data,
        positions=x,
        widths=0.6,
        showmeans=True,
        showmedians=True,
        showextrema=False,
    )
    for body in vp1['bodies']:
        body.set_facecolor(color)
        body.set_alpha(0.5)
        body.set_edgecolor("k")
        body.set_linewidth(0.5)
    vp1['cmeans'].set_color("k")
    vp1['cmeans'].set_linewidth(1)
    vp1['cmedians'].set_color("k")
    vp1['cmedians'].set_linewidth(1)
    vp1['cmedians'].set_linestyle("--")

    ax_elbo.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
    ax_elbo.set_xticks(x)
    ax_elbo.set_xticklabels([f"{cm}" for cm in cm_values], fontsize=10)
    ax_elbo.set_xlabel("Genome length (cM)", fontsize=11)
    ax_elbo.set_ylabel("ELBO(T_true) - ELBO(T_wrong)", fontsize=11)
    ax_elbo.set_title(f"{m_label}: delta-ELBO", fontsize=13, fontweight="bold")
    ax_elbo.grid(axis="y", alpha=0.15)

    # --- Bottom: w(T_true) ---
    vp2 = ax_w.violinplot(
        weight_data,
        positions=x,
        widths=0.6,
        showmeans=True,
        showmedians=True,
        showextrema=False,
    )
    for body in vp2['bodies']:
        body.set_facecolor(color)
        body.set_alpha(0.5)
        body.set_edgecolor("k")
        body.set_linewidth(0.5)
    vp2['cmeans'].set_color("k")
    vp2['cmeans'].set_linewidth(1)
    vp2['cmedians'].set_color("k")
    vp2['cmedians'].set_linewidth(1)
    vp2['cmedians'].set_linestyle("--")

    ax_w.axhline(0.5, color='k', ls='--', lw=1, alpha=0.5, label="Chance (0.50)")
    ax_w.set_xticks(x)
    ax_w.set_xticklabels([f"{cm}" for cm in cm_values], fontsize=10)
    ax_w.set_xlabel("Genome length (cM)", fontsize=11)
    ax_w.set_ylabel("w(T_true)", fontsize=11)
    ax_w.set_ylim(-0.05, 1.15)
    ax_w.set_title(f"{m_label}: w(T_true)", fontsize=13, fontweight="bold")
    ax_w.legend(fontsize=10)
    ax_w.grid(axis="y", alpha=0.15)

    # Per-cm misclassification rate below each boxplot
    for ci, cm_val in enumerate(cm_values):
        wd = weight_data[ci]
        if not np.all(np.isnan(wd)):
            pct = 100.0 * np.mean(wd < 0.5)
            ax_w.text(ci, -0.03, f"{pct:.0f}%", ha='center', va='top',
                      fontsize=8, color='red', fontweight='bold')

    # Label the bottom row
    ax_w.text(0.5, -0.12, "misclassified %", ha='center', va='top',
              fontsize=9, color='red', transform=ax_w.transAxes)

fig.suptitle(
    f"Topology swap (Nvarying): T_true = pure tree ((a,b),c),d) vs T_wrong = admixture on b\n"
    f"Ne(a)={NE_A}, Ne(rest)={NE_REST}",
    fontsize=13, fontweight="bold"
)

fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig("topology_swap_Nvarying.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: topology_swap_Nvarying.png")


# ====================================================================
# 6. Summary table
# ====================================================================
print("\n" + "=" * 90)
print("TOPOLOGY SWAP SUMMARY (Nvarying, all three models)")
print("=" * 90)
print(f"{'cm':>6s} | {'Model':>10s} | {'mean dELBO':>12s} | {'median':>10s} | "
      f"{'mean w(T1)':>10s} | {'% T1 wins':>10s} | {'n_ok':>5s}")
print("-" * 90)

for cm_val in cm_values:
    for m_label in model_labels:
        elbos = elbo_results[m_label][cm_val]
        if len(elbos) > 0:
            arr = np.array(elbos)
            diff = arr[:, 0] - arr[:, 1]
            max_e = np.max(arr, axis=1, keepdims=True)
            exp_e = np.exp(arr - max_e)
            w = exp_e / exp_e.sum(axis=1, keepdims=True)
            pct_win = 100.0 * np.mean(diff > 0)
            print(f"{cm_val:6d} | {m_label:>10s} | {np.mean(diff):+11.1f} | "
                  f"{np.median(diff):+9.1f} | "
                  f"{np.mean(w[:, 0]):10.3f} | {pct_win:9.1f}% | {len(elbos):5d}")
        else:
            print(f"{cm_val:6d} | {m_label:>10s} | {'N/A':>12s} | {'N/A':>10s} | "
                  f"{'N/A':>10s} | {'N/A':>10s} | {0:5d}")
