"""
Diagnostic: compare observed IBD fractions from simulation against
analytical predictions (int_p_L) evaluated at the TRUE parameter values.

If there's a systematic discrepancy, it explains the Ne bias.
"""

import numpy as np
from collections import defaultdict

from simulation_methods import simulate_msprime, build_ibd_stan_data
from ibd_jackknife import calculate_ibd_blocks_mrca, pool_multiple_simulations
from demography import DemographicTopology


# ====================================================================
# 0. Configuration (same as run_jackknife_analysis.py)
# ====================================================================
BLOCK_SIZE_CM = 50.0
CM_PER_UNIT = 1e-4
RECOMB_RATE = 1e-6
MUT_RATE = 1e-6
SAMPLES_PER_POP = {'a': 15, 'b': 15, 'c': 15, 'admix': 15}
N_SIMS = 5
SIM_CM_EACH = 500
SIM_SEQ_LEN = SIM_CM_EACH / CM_PER_UNIT

bins = [
    [0.5, 0.55], [0.55, 0.6], [0.6, 0.7], [0.7, 0.8], [0.8, 0.9], [0.9, 1.0],
    [1.0, 1.5], [1.5, 2.0], [2.0, 5.0], [5.0, 10.0], [10.0, 500.0]
]

# TRUE parameter values
TRUE_NE = 6000.0  # haploid effective N
TRUE_TIMES = [20, 40, 50, 200, 500]  # cumulative event times

# ====================================================================
# 1. Define topology (same as main scripts)
# ====================================================================
dem = DemographicTopology(['a', 'admix', 'b', 'c'])
dem.add_admixture_event('admix', 'admixPa1', 'admixPa2')
dem.add_merge_event('a', 'admixPa1', 'aAdmix')
dem.add_merge_event('b', 'admixPa2', 'bAdmix')
dem.add_merge_event('aAdmix', 'bAdmix', 'subAnc')
dem.add_merge_event('subAnc', 'c', 'anc')
for node in ['a','b','c','admix','admixPa1','admixPa2','aAdmix','bAdmix','subAnc','anc']:
    dem.set_node_ne(node, 3000)
dem.set_admixture_parameters('admix', 20, 0.75, 'admixPa1')
dem.set_merge_time('aAdmix', 40)
dem.set_merge_time('bAdmix', 50)
dem.set_merge_time('subAnc', 200)
dem.set_merge_time('anc', 500)
dem.finalize_root()


# ====================================================================
# 2. Analytical prediction: int_p_L in Python
# ====================================================================
def int_p_L(N, t1, t2, u, v):
    """Reproduce the Stan int_p_L function in Python."""
    if t2 - t1 < 1e-10:
        return 0.0

    cu = -(1.0/N + u/50.0)
    cv = -(1.0/N + v/50.0)

    eu1 = np.exp(-u * t1 / 50.0)
    eu2 = np.exp((t1 - t2) / N - u * t2 / 50.0)
    ev1 = np.exp(-v * t1 / 50.0)
    ev2 = np.exp((t1 - t2) / N - v * t2 / 50.0)

    cu2 = cu * cu
    cv2 = cv * cv

    k1 = u/(50.0*N) * (eu2 * (cu*t2 - 1.0)/cu2 - eu1 * (cu*t1 - 1.0)/cu2)
    k2 = v/(50.0*N) * (ev2 * (cv*t2 - 1.0)/cv2 - ev1 * (cv*t1 - 1.0)/cv2)

    k3 = 1.0/N * (eu2/cu - eu1/cu)
    k4 = 1.0/N * (ev2/cv - ev1/cv)

    return k1 - k2 + k3 - k4


def compute_predicted_ibd(dem, bins, N, T_max=100000):
    """
    Compute predicted IBD fractions using the same logic as the Stan model,
    evaluated at the TRUE parameter values.

    Returns dict: predicted[b_i][leaf_i, leaf_j]
    """
    matrices, events, admixture_map, admixture_map_id = \
        dem.get_topology_matrix_representation()

    n_nodes = len(dem.nodes)
    n_leaves = len(dem.initial_leaves)
    n_events = len(dem.ordered_events)

    # Build cumulative times from true event times
    cumulative_times = np.array(TRUE_TIMES, dtype=float)

    # Time intervals (durations)
    durations = np.diff(np.concatenate([[0], cumulative_times]))

    # Build parameter migration matrices (identity + admixture + merge)
    I = np.eye(n_nodes)

    # Start with identity for epoch 1 (before first event)
    param_mats = [I.copy()]  # index 0 = before first event

    # Build the admixture indices and fixed indices (matching Stan)
    admixture_indices_0based = list(admixture_map_id.keys())

    # Events are ordered: each event produces a migration matrix
    # The Stan code builds parameter_migration_matrices[1..n_events+1]
    # param_mats[0] = I (epoch 1, before any event)
    # param_mats[e] = migration matrix applied before epoch e+1

    for e_idx in range(n_events):
        if e_idx in admixture_map_id:
            # Admixture event
            info = admixture_map_id[e_idx]
            src = info['child']
            tgt1, tgt2 = info['parents']
            admix_mat = I.copy()
            admix_mat[src, tgt1] = 0.75  # TRUE admixture fraction
            admix_mat[src, tgt2] = 0.25
            admix_mat[src, src] = 0.0
            param_mats.append(admix_mat)
        else:
            param_mats.append(matrices[e_idx].copy())

    # Compute IBD fractions for each leaf pair
    # Leaf pairs: enumerate all (i, j) with i <= j
    predicted = {}
    for b_i in range(len(bins)):
        predicted[b_i] = np.zeros((n_leaves, n_leaves))

    for li in range(n_leaves):
        for lj in range(li, n_leaves):
            # Initialize W matrix
            W = np.zeros((n_nodes, n_nodes))
            W[li, lj] = 1.0

            for e in range(n_events + 1):
                # Migration event before epoch e (for e > 0)
                if e > 0:
                    M = param_mats[e]
                    W = M.T @ W @ M

                # Epoch boundaries
                if e == 0:
                    t_start = 0.0
                    t_end = cumulative_times[0]
                elif e == n_events:
                    t_start = cumulative_times[-1]
                    t_end = T_max
                else:
                    t_start = cumulative_times[e - 1]
                    t_end = cumulative_times[e]

                duration = t_end - t_start

                # IBD accumulation + survival update
                for a in range(n_nodes):
                    w_diag = W[a, a]
                    if w_diag > 1e-20:
                        for b_i in range(len(bins)):
                            u, v = bins[b_i]
                            predicted[b_i][li, lj] += w_diag * int_p_L(N, t_start, t_end, u, v)

                        # Survival: no coalescence in this epoch
                        W[a, a] = w_diag * np.exp(-duration / N)

            # Fill symmetric
            if li != lj:
                for b_i in range(len(bins)):
                    predicted[b_i][lj, li] = predicted[b_i][li, lj]

    return predicted


# ====================================================================
# 3. Simulate and compute observed IBD
# ====================================================================
print("="*60)
print("Step 1: Simulating and computing observed IBD")
print("="*60)

packed_list = []
n_blocks_list = []
pair_info = None

for sim_i in range(N_SIMS):
    print(f"\n--- Simulation {sim_i + 1}/{N_SIMS} ({SIM_CM_EACH} cM) ---")
    mts = simulate_msprime(
        dem,
        sequence_length=SIM_SEQ_LEN,
        recombination_rate=RECOMB_RATE,
        mutation_rate=MUT_RATE,
        samples_per_pop=SAMPLES_PER_POP,
        seed=42 + sim_i,
    )

    packed, n_blocks, pop_samples, pop_ids, pi = calculate_ibd_blocks_mrca(
        mts,
        bins=bins,
        block_size_cm=BLOCK_SIZE_CM,
        cm_per_unit=CM_PER_UNIT,
    )
    packed_list.append(packed)
    n_blocks_list.append(n_blocks)
    if pair_info is None:
        pair_info = pi

pooled, total_blocks = pool_multiple_simulations(packed_list, n_blocks_list)

# Compute observed mean over ALL blocks (no resampling, using full data)
num_pops = len(pop_ids)
n_bins = len(bins)
observed = {}

for b_i in range(n_bins):
    obs_mat = np.zeros((num_pops, num_pops))
    for i in range(num_pops):
        for j in range(i, num_pops):
            arr = pooled[(b_i, i, j)]
            if arr.shape[0] > 0:
                pair_fracs = arr.sum(axis=1) / total_blocks
                obs_mat[i, j] = obs_mat[j, i] = pair_fracs.mean()
    observed[b_i] = obs_mat

# Print number of haplotypes per population
print("\n" + "="*60)
print("Haplotype counts per population:")
for pid in pop_ids:
    print(f"  Pop {pid}: {len(pop_samples[pid])} haplotypes")
print("="*60)

# ====================================================================
# 4. Compute predicted IBD at true parameters
# ====================================================================
print("\nStep 2: Computing analytical predictions at true parameters")
predicted = compute_predicted_ibd(dem, bins, TRUE_NE, T_max=100000)

# ====================================================================
# 5. Compare observed vs predicted
# ====================================================================
print("\n" + "="*70)
print(f"{'Bin':>12s} | {'Pop pair':>10s} | {'Observed':>12s} | {'Predicted':>12s} | {'Ratio':>8s} | {'Diff%':>8s}")
print("-"*70)

all_ratios = []
weighted_ratios = []

for b_i in range(n_bins):
    bin_str = f"[{bins[b_i][0]:.2f},{bins[b_i][1]:.0f}]" if bins[b_i][1] > 10 else f"[{bins[b_i][0]:.2f},{bins[b_i][1]:.2f}]"
    for i in range(num_pops):
        for j in range(i, num_pops):
            obs = observed[b_i][i, j]
            pred = predicted[b_i][i, j]

            if pred > 1e-10:
                ratio = obs / pred
                diff_pct = (obs - pred) / pred * 100
                all_ratios.append(ratio)
                weighted_ratios.append((ratio, pred))
            else:
                ratio = float('nan')
                diff_pct = float('nan')

            if obs > 1e-8 or pred > 1e-8:
                print(f"{bin_str:>12s} | ({i},{j}){' '*4} | {obs:12.6e} | {pred:12.6e} | {ratio:8.4f} | {diff_pct:+7.2f}%")

# Summary statistics
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
valid_ratios = [r for r in all_ratios if not np.isnan(r)]
print(f"Mean obs/pred ratio:   {np.mean(valid_ratios):.4f}")
print(f"Median obs/pred ratio: {np.median(valid_ratios):.4f}")

# Weighted mean (by predicted value)
weights = np.array([w for _, w in weighted_ratios])
ratios = np.array([r for r, _ in weighted_ratios])
valid = ~np.isnan(ratios)
wmean = np.average(ratios[valid], weights=weights[valid])
print(f"Weighted mean ratio:   {wmean:.4f}")
print(f"\nIf ratio > 1: observed > predicted (model underestimates IBD → would bias Ne downward)")
print(f"If ratio < 1: observed < predicted (model overestimates IBD → would bias Ne upward)")
