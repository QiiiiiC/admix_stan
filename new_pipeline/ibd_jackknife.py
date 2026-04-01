"""
Block-bootstrap IBD resampling with delete-one-haplotype jackknife variance.

Same block-bootstrap strategy as bootstrap_ibd.py for genome resampling,
but uses a jackknife over haplotypes (the independent sampling units) instead
of a pairs bootstrap for variance estimation. This correctly accounts for
the positive correlation between pairs sharing a haplotype.

The mean over pairs is a U-statistic of degree 2. The delete-one jackknife
captures the ~4x variance inflation from the U-statistic structure:
    Var(mean over pairs) ≈ 4 * Var(g_u) / n
where g_u = average IBD of haplotype u with all partners.
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
    Pairs are packed in deterministic order with known haplotype identity.

    Returns
    -------
    packed : dict
        packed[(b_i, i, j)] = np.array of shape (num_pairs, n_blocks).
        Rows are in deterministic order matching pair_info.

    n_blocks : int

    pop_samples : dict  {pop_id: [sample_node_ids]}

    pop_ids : np.array

    pair_info : dict
        pair_info[(i, j)] = dict with keys:
            'pair_a': np.array of local haplotype indices (first member)
            'pair_b': np.array of local haplotype indices (second member)
            'n_haps_i': number of haplotypes in pop i
            'n_haps_j': number of haplotypes in pop j
            'is_within': bool, True if i == j
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

    # --- Phase 1: collect into sparse dicts (same as bootstrap_ibd.py) ---
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

    # --- Phase 2: pack with deterministic pair ordering + pair identity ---
    packed = {}
    pair_info = {}

    for i in range(num_pops):
        for j in range(i, num_pops):
            pi, pj = pop_ids[i], pop_ids[j]
            key = (min(pi, pj), max(pi, pj))

            haps_i = pop_samples[pop_ids[i]]
            haps_j = pop_samples[pop_ids[j]]
            n_haps_i = len(haps_i)
            n_haps_j = len(haps_j)

            # Enumerate all pairs in deterministic order
            all_pair_a = []
            all_pair_b = []
            all_pair_keys = []

            if i == j:
                for a in range(n_haps_i):
                    for b in range(a + 1, n_haps_i):
                        all_pair_a.append(a)
                        all_pair_b.append(b)
                        all_pair_keys.append(tuple(sorted((haps_i[a], haps_i[b]))))
            else:
                for a in range(n_haps_i):
                    for b in range(n_haps_j):
                        all_pair_a.append(a)
                        all_pair_b.append(b)
                        all_pair_keys.append(tuple(sorted((haps_i[a], haps_j[b]))))

            num_pairs = len(all_pair_keys)

            pair_info[(i, j)] = {
                'pair_a': np.array(all_pair_a, dtype=np.int32),
                'pair_b': np.array(all_pair_b, dtype=np.int32),
                'n_haps_i': n_haps_i,
                'n_haps_j': n_haps_j,
                'is_within': i == j,
            }

            # Fill packed arrays for each bin
            for b_i in range(len(bins)):
                observed_dict = raw[b_i].get(key, {})
                arr = np.zeros((num_pairs, n_blocks))

                for row, pk in enumerate(all_pair_keys):
                    if pk in observed_dict:
                        arr[row] = observed_dict[pk]

                packed[(b_i, i, j)] = arr

    print("Packing complete (deterministic pair ordering).")
    return packed, n_blocks, dict(pop_samples), pop_ids, pair_info


def pool_multiple_simulations(packed_list, n_blocks_list):
    """
    Concatenate block arrays from multiple independent simulations.
    """
    total_blocks = sum(n_blocks_list)
    keys = list(packed_list[0].keys())

    pooled = {}
    for key in keys:
        arrays = [p[key] for p in packed_list]
        pooled[key] = np.hstack(arrays)

    print(f"Pooled {len(packed_list)} simulations: {total_blocks} total blocks")
    return pooled, total_blocks


def _jackknife_variance_within(pair_fracs, pair_a, pair_b, n_haps):
    """
    Delete-one-haplotype jackknife for within-population pair mean.

    For n haplotypes and P = n(n-1)/2 pairs, the mean is a one-sample
    U-statistic. The jackknife pseudovalue for haplotype h is:
        theta_h = n * mu - (n-1) * mu_(-h)
    and Var(mu) = sum(theta_h - mu)^2 / (n*(n-1)).

    This equals approximately 4 * Var(g_h) / n, where g_h is the mean
    IBD of haplotype h with all partners.
    """
    n = n_haps
    P = len(pair_fracs)
    S = pair_fracs.sum()
    mu = S / P

    # For each haplotype h, sum of pair_fracs for pairs involving h
    hap_sums = np.zeros(n)
    np.add.at(hap_sums, pair_a, pair_fracs)
    np.add.at(hap_sums, pair_b, pair_fracs)

    # Each haplotype appears in exactly (n-1) within-pop pairs
    count_per_hap = n - 1
    P_remaining = P - count_per_hap

    # Leave-one-out means
    leave_out_means = (S - hap_sums) / P_remaining

    # Pseudovalues
    thetas = n * mu - (n - 1) * leave_out_means

    # Jackknife variance (mean(thetas) = mu exactly)
    var = np.sum((thetas - mu) ** 2) / (n * (n - 1))
    return var


def _jackknife_variance_between(pair_fracs, pair_a, pair_b, n_haps_i, n_haps_j):
    """
    Delete-one-haplotype jackknife for between-population pair mean.

    For n_i x n_j pairs (two-sample U-statistic), jackknife from both
    sides independently and sum the variance contributions:
        Var(mu) = Var_i + Var_j
    """
    P = len(pair_fracs)
    S = pair_fracs.sum()
    mu = S / P

    # --- Jackknife from pop_i side ---
    hap_sums_i = np.zeros(n_haps_i)
    np.add.at(hap_sums_i, pair_a, pair_fracs)

    # Each haplotype in pop_i appears in n_j pairs
    P_remaining_i = P - n_haps_j
    leave_out_means_i = (S - hap_sums_i) / P_remaining_i
    thetas_i = n_haps_i * mu - (n_haps_i - 1) * leave_out_means_i
    var_i = np.sum((thetas_i - mu) ** 2) / (n_haps_i * (n_haps_i - 1))

    # --- Jackknife from pop_j side ---
    hap_sums_j = np.zeros(n_haps_j)
    np.add.at(hap_sums_j, pair_b, pair_fracs)

    # Each haplotype in pop_j appears in n_i pairs
    P_remaining_j = P - n_haps_i
    leave_out_means_j = (S - hap_sums_j) / P_remaining_j
    thetas_j = n_haps_j * mu - (n_haps_j - 1) * leave_out_means_j
    var_j = np.sum((thetas_j - mu) ** 2) / (n_haps_j * (n_haps_j - 1))

    return var_i + var_j


def resample_ibd_with_jackknife_variance(
    packed: dict,
    n_blocks_total: int,
    pop_ids: np.ndarray,
    pair_info: dict,
    bins: List[Tuple[float, float]],
    target_cm: float,
    block_size_cm: float = 50.0,
    rng: np.random.Generator = None,
    chosen_blocks: np.ndarray = None,
):
    """
    Block-bootstrap one IBD fraction matrix with jackknife variance estimation.

    Block bootstrap draws m = target_cm / block_size_cm blocks with replacement
    to simulate a genome of length target_cm. The jackknife over haplotypes
    gives the variance of the mean IBD fraction across pairs.

    If chosen_blocks is provided, those block indices are used directly
    (for synchronized resampling with SNP data).
    """
    if rng is None:
        rng = np.random.default_rng()

    if chosen_blocks is None:
        n_blocks_needed = max(1, int(round(target_cm / block_size_cm)))
        chosen_blocks = rng.integers(0, n_blocks_total, size=n_blocks_needed)

    n_blocks_needed = len(chosen_blocks)

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

                info = pair_info[(i, j)]

                if info['is_within']:
                    if info['n_haps_i'] > 2:
                        var_val = _jackknife_variance_within(
                            pair_fracs, info['pair_a'], info['pair_b'],
                            info['n_haps_i'],
                        )
                    else:
                        var_val = se_floor ** 2
                else:
                    if info['n_haps_i'] > 1 and info['n_haps_j'] > 1:
                        var_val = _jackknife_variance_between(
                            pair_fracs, info['pair_a'], info['pair_b'],
                            info['n_haps_i'], info['n_haps_j'],
                        )
                    else:
                        var_val = se_floor ** 2

                mean_mat[i, j] = mean_mat[j, i] = mean_val
                var_mat[i, j] = var_mat[j, i] = max(var_val, se_floor ** 2)

        ibd_mean[b_i] = mean_mat
        ibd_var[b_i] = var_mat

    return ibd_mean, ibd_var
