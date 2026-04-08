"""
Topology swap experiment: Can the mixed model distinguish the temporal order
of post-admixture merges?

Data is simulated under T1 (true):
  admix = 75% from a's lineage + 25% from b's lineage
  Event order (backward in time):
    1. admixture: admix → admixPa1 (75%) + admixPa2 (25%)
    2. a + admixPa1 (major) → aAdmix   (first merge, time 40)
    3. b + admixPa2 (minor) → bAdmix   (second merge, time 50)
    4. aAdmix + bAdmix → subAnc         (time 1000)
    5. subAnc + c → anc                  (time 1500)

T_rev (reversed merge order):
  Same populations merge with same parents, but b merges first:
    1. admixture: admix → admixPa1 (75%) + admixPa2 (25%)
    2. b + admixPa2 (minor) → bAdmix   (first merge)
    3. a + admixPa1 (major) → aAdmix   (second merge)
    4. bAdmix + aAdmix → subAnc
    5. subAnc + c → anc

Only the mixed model is used.
Plot: ELBO(T1) - ELBO(T_rev) vs genome length.
"""

import sys, os
_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _PARENT)
_MODELS = os.path.join(_PARENT, "models")

import numpy as np
import matplotlib.pyplot as plt
import re
from contextlib import redirect_stdout

from simulation_methods import (
    simulate_msprime,
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
SAMPLES_PER_POP = {'a': 15, 'b': 15, 'c': 15, 'admix': 15}

N_SIMS = 50
SIM_CM_EACH = 50
SIM_SEQ_LEN = SIM_CM_EACH / CM_PER_UNIT

cm_values = [50, 100, 150, 200, 250, 300, 400, 500]

bins = [
    [0.5, 0.55], [0.55, 0.6], [0.6, 0.65], [0.65, 0.7], [0.7, 0.8], [0.8, 0.9], [0.9, 1.0],
    [1.0, 1.5], [1.5, 2.0], [2.0, 5.0], [5.0, 8.0], [8.0, 20.0], [20.0, BLOCK_SIZE_CM]
]

SNP_CUTOFF_TIME = 1500
SNP_MIN_MAF = 0.05


# ====================================================================
# 1. Define topologies
# ====================================================================
print("=" * 60)
print("Defining topologies")
print("=" * 60)

# --- T1 (true): a merges with major parent FIRST, b with minor SECOND ---
dem_sim = DemographicTopology(['a', 'admix', 'b', 'c'])
dem_sim.add_admixture_event('admix', 'admixPa1', 'admixPa2')
dem_sim.add_merge_event('a', 'admixPa1', 'aAdmix')       # first merge (time 40)
dem_sim.add_merge_event('b', 'admixPa2', 'bAdmix')       # second merge (time 50)
dem_sim.add_merge_event('aAdmix', 'bAdmix', 'subAnc')
dem_sim.add_merge_event('subAnc', 'c', 'anc')
for node in dem_sim.nodes:
    dem_sim.set_node_ne(node, 3000)
dem_sim.set_admixture_parameters('admix', 20, 0.75, 'admixPa1')
dem_sim.set_merge_time('aAdmix', 40)
dem_sim.set_merge_time('bAdmix', 50)
dem_sim.set_merge_time('subAnc', 1000)
dem_sim.set_merge_time('anc', 1500)
dem_sim.finalize_root()

dem_T1 = dem_sim

# --- T_rev: b merges with minor parent FIRST, a with major SECOND ---
dem_rev = DemographicTopology(['a', 'admix', 'b', 'c'])
dem_rev.add_admixture_event('admix', 'admixPa1', 'admixPa2')
dem_rev.add_merge_event('b', 'admixPa2', 'bAdmix')       # first merge (b+minor)
dem_rev.add_merge_event('a', 'admixPa1', 'aAdmix')       # second merge (a+major)
dem_rev.add_merge_event('bAdmix', 'aAdmix', 'subAnc')
dem_rev.add_merge_event('subAnc', 'c', 'anc')
dem_rev.finalize_root()

topology_dems = [dem_T1, dem_rev]
topology_labels = ["T1 (true: a+major first)", "T_rev (b+minor first)"]

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
# 3. Compile Stan model (mixed only)
# ====================================================================
mixed_stan = CmdStanModel(stan_file=os.path.join(_MODELS, "mixed_model_Nfixed.stan"))


# ====================================================================
# 4. Run replicates: fit T1 vs T_rev with mixed model
# ====================================================================
elbo_all = {cm: [] for cm in cm_values}

rng = np.random.default_rng(seed=2025)
_devnull = open(os.devnull, 'w')

for cm_val in cm_values:
    print(f"\n{'='*60}")
    print(f"  cm = {cm_val} cM  ({N_REPLICATES} reps × 2 topologies)")
    print(f"{'='*60}")

    n_blocks_needed = max(1, int(round(cm_val / BLOCK_SIZE_CM)))

    for rep in range(N_REPLICATES):
        if (rep + 1) % 10 == 0:
            print(f"  replicate {rep + 1}/{N_REPLICATES}")

        chosen_blocks = rng.choice(total_blocks, size=n_blocks_needed, replace=False)

        # --- Resample data (shared across topologies) ---
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

        # --- Fit each topology ---
        rep_elbos = [np.nan, np.nan]

        for ti, dem_infer in enumerate(topology_dems):
            n_events = len(dem_infer.ordered_events)
            n_admix = dem_infer.n_admix

            init = {
                "times": [100.0] * n_events,
                "admixture_fractions": [0.5] * n_admix,
                "effective_N": 10000.0,
                "kappa_snp": 1.0,
            }

            try:
                with redirect_stdout(_devnull):
                    sd = build_mixed_stan_data(
                        dem_infer, ibd_mean, ibd_var, bins, w_hat, w_se,
                        T_max=100000, cm=cm_val,
                    )
                fit = mixed_stan.pathfinder(
                    data=sd, inits=init, show_console=False, psis_resample=True,
                )
                rep_elbos[ti] = extract_elbo(fit)
            except Exception as e:
                if rep < 3:
                    print(f"    [WARN] {topology_labels[ti]} rep {rep+1}: {e}")

        if all(np.isfinite(rep_elbos)):
            elbo_all[cm_val].append(tuple(rep_elbos))

_devnull.close()


# ====================================================================
# 5. Plot: two panels — ΔELBO and w(T1)
# ====================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))

diff_data = []
weight_data = []
for cm_val in cm_values:
    elbos = elbo_all[cm_val]
    if len(elbos) > 0:
        arr = np.array(elbos)  # (n_reps, 2)
        diff_data.append(arr[:, 0] - arr[:, 1])
        # softmax weight for T1: exp(e1) / (exp(e1) + exp(e2))
        max_e = np.max(arr, axis=1, keepdims=True)
        exp_e = np.exp(arr - max_e)
        w = exp_e / exp_e.sum(axis=1, keepdims=True)
        weight_data.append(w[:, 0])
    else:
        diff_data.append(np.array([np.nan]))
        weight_data.append(np.array([np.nan]))

x = np.arange(len(cm_values))

# --- Panel 1: ΔELBO ---
bp1 = ax1.boxplot(
    diff_data,
    positions=x,
    widths=0.5,
    patch_artist=True,
    showfliers=False,
    medianprops=dict(color="k", linewidth=1.5),
    whiskerprops=dict(color="#E8A838", linewidth=1),
    capprops=dict(color="#E8A838", linewidth=1),
)
for patch in bp1['boxes']:
    patch.set_facecolor("#E8A838")
    patch.set_alpha(0.4)

means = [np.nanmean(d) for d in diff_data]
ax1.plot(x, means, color="#E8A838", linewidth=2, marker='o', markersize=5, zorder=5,
         label="Mean ΔELBO")

ax1.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
ax1.set_xticks(x)
ax1.set_xticklabels([f"{cm}" for cm in cm_values])
ax1.set_xlabel("Genome length (cM)")
ax1.set_ylabel("ΔELBO  (T_true − T_rev)")
ax1.set_title("ELBO difference\nPositive = true order preferred")
ax1.legend(fontsize=9)
ax1.grid(axis="y", alpha=0.15)

# --- Panel 2: w(T1) ---
bp2 = ax2.boxplot(
    weight_data,
    positions=x,
    widths=0.5,
    patch_artist=True,
    showfliers=True,
    flierprops=dict(marker='o', markersize=3, markerfacecolor="#E8A838", alpha=0.4),
    medianprops=dict(color="k", linewidth=1.5),
    whiskerprops=dict(color="#E8A838", linewidth=1),
    capprops=dict(color="#E8A838", linewidth=1),
)
for patch in bp2['boxes']:
    patch.set_facecolor("#E8A838")
    patch.set_alpha(0.4)

mean_w = [np.nanmean(d) for d in weight_data]
ax2.plot(x, mean_w, color="#E8A838", linewidth=2, marker='o', markersize=5, zorder=5,
         label="Mean w(T1)")

ax2.axhline(0.5, color='k', ls='--', lw=1, alpha=0.5, label="Chance (0.50)")
ax2.set_xticks(x)
ax2.set_xticklabels([f"{cm}" for cm in cm_values])
ax2.set_xlabel("Genome length (cM)")
ax2.set_ylabel("w(T_true)")
ax2.set_title("Model weight for true topology\nw = exp(ELBO) / Σ exp(ELBO)")
ax2.set_ylim(-0.05, 1.05)
ax2.legend(fontsize=9)
ax2.grid(axis="y", alpha=0.15)

# Compute overall misclassification rate (w(T1) < 0.5 across all cm_values)
all_weights = np.concatenate([d for d in weight_data if not np.all(np.isnan(d))])
misclass_pct = 100.0 * np.mean(all_weights < 0.5)

fig.suptitle("Mixed model: can we distinguish merge order?\n"
             "True: a+major first, then b+minor  |  Rev: b+minor first, then a+major\n"
             f"Misclassified (w(T_true) < 0.5): {misclass_pct:.1f}%",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig("topology_swap_elbo.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: topology_swap_elbo.png")


# ====================================================================
# 6. Summary table
# ====================================================================
print("\n" + "=" * 70)
print("TOPOLOGY SWAP SUMMARY (Mixed model)")
print("=" * 70)
print(f"{'cm':>6s} | {'mean ΔELBO':>12s} | {'median':>10s} | {'mean w(T1)':>10s} | "
      f"{'% T1 wins':>10s} | {'n_ok':>5s}")
print("-" * 70)

for cm_val in cm_values:
    elbos = elbo_all[cm_val]
    if len(elbos) > 0:
        arr = np.array(elbos)
        diff = arr[:, 0] - arr[:, 1]
        max_e = np.max(arr, axis=1, keepdims=True)
        exp_e = np.exp(arr - max_e)
        w = exp_e / exp_e.sum(axis=1, keepdims=True)
        pct_win = 100.0 * np.mean(diff > 0)
        print(f"{cm_val:6d} | {np.mean(diff):+11.1f} | {np.median(diff):+9.1f} | "
              f"{np.mean(w[:, 0]):10.3f} | {pct_win:9.1f}% | {len(elbos):5d}")
    else:
        print(f"{cm_val:6d} | {'N/A':>12s} | {'N/A':>10s} | {'N/A':>10s} | "
              f"{'N/A':>10s} | {0:5d}")
