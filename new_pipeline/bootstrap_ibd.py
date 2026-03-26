"""
Block-bootstrap IBD resampling module (vectorized).

Strategy:
  1. Simulate ONCE with a large genome.
  2. Compute per-pair, per-block IBD contributions (expensive, done once).
  3. Pack everything into dense arrays so resampling is pure numpy indexing.
  4. Each resample: pick block indices, sum columns, take mean. No Python loops over pairs.
"""

import numpy as np
from collections import defaultdict
from typing import List, Tuple


def calculate_ibd_blocks_mrca(
    ts,
    bins: List[Tuple[float, float]],
    block_size_cm: float = 50.0,
    cm_per_unit: float = 1e-6,
):
    """
    Compute per-block, per-pair IBD contributions using the MRCA-span method.
    Segments are binned by FULL length, contributions distributed across blocks.

    Returns
    -------
    packed : dict
        packed[(b_i, i, j)] = np.array(num_pairs, n_blocks)
        Dense matrix including zero-IBD pairs. i <= j are pop indices.

    n_blocks : int

    pop_samples : dict  {pop_id: [sample_node_ids]}

    pop_ids : np.array
    """
    sample_nodes = ts.samples()
    num_samples = len(sample_nodes)
    node_to_pop = ts.nodes_population[sample_nodes]
    pop_ids = np.unique(node_to_pop)
    num_pops = len(pop_ids)

    pop_samples = defaultdict(list)
    for u in sample_nodes:
        pop_samples[node_to_pop[u]].append(u)

    genome_length_bp = ts.sequence_length
    genome_length_cm = genome_length_bp * cm_per_unit
    block_size_bp = block_size_cm / cm_per_unit
    n_blocks = int(np.floor(genome_length_cm / block_size_cm))

    if n_blocks == 0:
        raise ValueError(
            f"Genome ({genome_length_cm:.1f} cM) shorter than one block ({block_size_cm} cM)."
        )

    used_length_bp = n_blocks * block_size_bp
    block_boundaries_bp = [i * block_size_bp for i in range(n_blocks + 1)]

    print(f"Genome: {genome_length_cm:.1f} cM -> {n_blocks} blocks of {block_size_cm} cM "
          f"(using {n_blocks * block_size_cm:.1f} cM)")

    # --- Phase 1: collect into sparse dicts (same as before) ---
    raw = {
        b_i: defaultdict(lambda: defaultdict(lambda: np.zeros(n_blocks)))
        for b_i in range(len(bins))
    }

    current_state = {}
    pairs = []

    tree_iter = ts.trees()
    first_tree = next(tree_iter)

    for i in range(num_samples):
        for j in range(i + 1, num_samples):
            u, v = sample_nodes[i], sample_nodes[j]
            pairs.append((u, v))
            mrca = first_tree.mrca(u, v)
            current_state[(u, v)] = (mrca, first_tree.interval.left)

    print(f"Scanning trees for MRCA segments ({len(pairs)} pairs)...")

    def assign_segment_to_blocks(u, v, seg_left_bp, seg_right_bp):
        seg_left_bp = max(seg_left_bp, 0.0)
        seg_right_bp = min(seg_right_bp, used_length_bp)
        if seg_right_bp <= seg_left_bp:
            return

        full_seg_cm = (seg_right_bp - seg_left_bp) * cm_per_unit

        target_bin = -1
        for b_i, (min_len, max_len) in enumerate(bins):
            if min_len < full_seg_cm <= max_len:
                target_bin = b_i
                break
        if target_bin == -1:
            return

        p_u = node_to_pop[u]
        p_v = node_to_pop[v]
        p_i, p_j = sorted((p_u, p_v))
        pair_key = tuple(sorted((u, v)))

        first_block = int(seg_left_bp // block_size_bp)
        last_block = min(int(seg_right_bp // block_size_bp), n_blocks - 1)
        if seg_right_bp <= block_boundaries_bp[last_block] and last_block > first_block:
            last_block -= 1

        for blk in range(first_block, last_block + 1):
            blk_left = block_boundaries_bp[blk]
            blk_right = block_boundaries_bp[blk + 1]
            overlap_left = max(seg_left_bp, blk_left)
            overlap_right = min(seg_right_bp, blk_right)
            overlap_cm = (overlap_right - overlap_left) * cm_per_unit
            if overlap_cm <= 0:
                continue
            frac = overlap_cm / block_size_cm
            raw[target_bin][(p_i, p_j)][pair_key][blk] += frac

    for tree in tree_iter:
        current_left = tree.interval.left
        for (u, v) in pairs:
            new_mrca = tree.mrca(u, v)
            old_mrca, start_pos = current_state[(u, v)]
            if new_mrca != old_mrca:
                assign_segment_to_blocks(u, v, start_pos, current_left)
                current_state[(u, v)] = (new_mrca, current_left)

    for (u, v) in pairs:
        old_mrca, start_pos = current_state[(u, v)]
        assign_segment_to_blocks(u, v, start_pos, ts.sequence_length)

    print("Block-level IBD computation complete. Packing into dense arrays...")

    # --- Phase 2: pack into dense arrays ---
    # For each (bin, pop_i, pop_j), create array of shape (num_pairs, n_blocks)
    # including zero-IBD pairs as zero rows.
    packed = {}

    for b_i in range(len(bins)):
        for i in range(num_pops):
            for j in range(i, num_pops):
                pi, pj = pop_ids[i], pop_ids[j]
                key = (min(pi, pj), max(pi, pj))

                if i == j:
                    n_i = len(pop_samples[pop_ids[i]])
                    num_pairs = n_i * (n_i - 1) // 2
                else:
                    num_pairs = len(pop_samples[pop_ids[i]]) * len(pop_samples[pop_ids[j]])

                if num_pairs == 0:
                    packed[(b_i, i, j)] = np.zeros((0, n_blocks))
                    continue

                observed_dict = raw[b_i].get(key, {})
                n_observed = len(observed_dict)
                n_zeros = num_pairs - n_observed

                if n_observed > 0:
                    obs_arr = np.array(list(observed_dict.values()))  # (n_observed, n_blocks)
                    if n_zeros > 0:
                        arr = np.vstack([obs_arr, np.zeros((n_zeros, n_blocks))])
                    else:
                        arr = obs_arr
                else:
                    arr = np.zeros((num_pairs, n_blocks))

                packed[(b_i, i, j)] = arr

    print("Packing complete.")

    return packed, n_blocks, dict(pop_samples), pop_ids


def resample_ibd_with_bootstrap_variance(
    packed: dict,
    n_blocks_total: int,
    pop_samples: dict,
    pop_ids: np.ndarray,
    bins: List[Tuple[float, float]],
    target_cm: float,
    block_size_cm: float = 50.0,
    n_var_bootstraps: int = 200,
    rng: np.random.Generator = None,
):
    """
    Block-bootstrap one IBD fraction matrix. All operations are vectorized numpy.

    Complexity: O(n_bins * n_pop_pairs) array ops, no Python loops over sample pairs.
    """
    if rng is None:
        rng = np.random.default_rng()

    n_blocks_needed = max(1, int(round(target_cm / block_size_cm)))
    chosen_blocks = rng.integers(0, n_blocks_total, size=n_blocks_needed)

    num_pops = len(pop_ids)
    n_bins = len(bins)
    se_floor = 1e-8

    ibd_mean = {}
    ibd_var = {}

    for b_i in range(n_bins):
        mean_mat = np.zeros((num_pops, num_pops))
        var_mat = np.zeros((num_pops, num_pops))

        for i in range(num_pops):
            for j in range(i, num_pops):
                arr = packed[(b_i, i, j)]  # (num_pairs, n_blocks)
                num_pairs = arr.shape[0]

                if num_pairs == 0:
                    continue

                # Vectorized: select columns, sum across chosen blocks, get per-pair mean
                # arr[:, chosen_blocks] is (num_pairs, n_blocks_needed)
                pair_fracs = arr[:, chosen_blocks].sum(axis=1) / n_blocks_needed  # (num_pairs,)

                mean_val = pair_fracs.mean()

                # Bootstrap variance over pairs
                if n_var_bootstraps > 0 and num_pairs > 1:
                    # Draw bootstrap indices: (n_var_bootstraps, num_pairs)
                    boot_idx = rng.integers(0, num_pairs, size=(n_var_bootstraps, num_pairs))
                    boot_means = pair_fracs[boot_idx].mean(axis=1)  # (n_var_bootstraps,)
                    var_val = boot_means.var()
                else:
                    var_val = se_floor ** 2

                mean_mat[i, j] = mean_mat[j, i] = mean_val
                var_mat[i, j] = var_mat[j, i] = max(var_val, se_floor ** 2)

        ibd_mean[b_i] = mean_mat
        ibd_var[b_i] = var_mat

    return ibd_mean, ibd_var


def resample_ibd_matrices(
    packed: dict,
    n_blocks_total: int,
    pop_samples: dict,
    pop_ids: np.ndarray,
    bins: List[Tuple[float, float]],
    target_cm: float,
    block_size_cm: float = 50.0,
    rng: np.random.Generator = None,
):
    """
    Fast resample without variance estimation (floor variance).
    """
    if rng is None:
        rng = np.random.default_rng()

    n_blocks_needed = max(1, int(round(target_cm / block_size_cm)))
    chosen_blocks = rng.integers(0, n_blocks_total, size=n_blocks_needed)

    num_pops = len(pop_ids)
    n_bins = len(bins)

    ibd_mean = {}
    ibd_var = {}

    for b_i in range(n_bins):
        mean_mat = np.zeros((num_pops, num_pops))
        var_mat = np.full((num_pops, num_pops), 1e-12)

        for i in range(num_pops):
            for j in range(i, num_pops):
                arr = packed[(b_i, i, j)]
                num_pairs = arr.shape[0]
                if num_pairs == 0:
                    continue

                pair_fracs = arr[:, chosen_blocks].sum(axis=1) / n_blocks_needed
                mean_mat[i, j] = mean_mat[j, i] = pair_fracs.mean()

        ibd_mean[b_i] = mean_mat
        ibd_var[b_i] = var_mat

    return ibd_mean, ibd_var

