"""
Block-bootstrap IBD resampling module (vectorized).

Strategy:
  1. Simulate MULTIPLE independent genomes (e.g. 5 x 500 cM).
  2. Compute per-pair, per-block IBD contributions for each (expensive, done once).
  3. Pool all blocks across simulations into one array.
  4. Resample blocks with replacement — drawing from multiple realizations
     avoids converging to a single simulation's idiosyncratic coalescent history.
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

    # --- Phase 1: collect into sparse dicts ---
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
                    obs_arr = np.array(list(observed_dict.values()))
                    if n_zeros > 0:
                        arr = np.vstack([obs_arr, np.zeros((n_zeros, n_blocks))])
                    else:
                        arr = obs_arr
                else:
                    arr = np.zeros((num_pairs, n_blocks))

                packed[(b_i, i, j)] = arr

    print("Packing complete.")

    return packed, n_blocks, dict(pop_samples), pop_ids


def pool_multiple_simulations(packed_list, n_blocks_list):
    """
    Concatenate block arrays from multiple independent simulations.

    Parameters
    ----------
    packed_list : list of dict
        Each element is the `packed` output from calculate_ibd_blocks_mrca.
        All must have the same keys (same topology, same bins, same sample sizes).

    n_blocks_list : list of int
        Number of blocks from each simulation.

    Returns
    -------
    pooled : dict
        Same structure as packed, but with columns concatenated across simulations.

    total_blocks : int
        Total number of blocks across all simulations.
    """
    total_blocks = sum(n_blocks_list)
    keys = list(packed_list[0].keys())

    pooled = {}
    for key in keys:
        arrays = [p[key] for p in packed_list]
        # All should have same number of rows (same num_pairs)
        # Concatenate along columns (blocks axis)
        pooled[key] = np.hstack(arrays)

    print(f"Pooled {len(packed_list)} simulations: {total_blocks} total blocks")
    return pooled, total_blocks


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
                arr = packed[(b_i, i, j)]
                num_pairs = arr.shape[0]

                if num_pairs == 0:
                    continue

                pair_fracs = arr[:, chosen_blocks].sum(axis=1) / n_blocks_needed

                mean_val = pair_fracs.mean()

                if n_var_bootstraps > 0 and num_pairs > 1:
                    boot_idx = rng.integers(0, num_pairs, size=(n_var_bootstraps, num_pairs))
                    boot_means = pair_fracs[boot_idx].mean(axis=1)
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

