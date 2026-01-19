import numpy as np
import tskit
from itertools import combinations_with_replacement

def calculate_treemix_covariance(allele_freq_matrix, block_size_snps=500, min_maf=0.0):
    """
    Calculates the TreeMix covariance matrix and standard error.
    
    Parameters:
    -----------
    allele_freq_matrix : np.ndarray
        Dimension (n_snps x n_pops). Values are allele frequencies (0.0 to 1.0).
    block_size_snps : int
        Number of SNPs per block for jackknife standard error estimation.
    min_maf : float
        Minimum Minor Allele Frequency to retain a SNP (default 0.0 to keep all valid).
        
    Returns:
    --------
    W_observed : np.ndarray (n_pops x n_pops)
        The observed covariance matrix.
    W_se : np.ndarray (n_pops x n_pops)
        The standard error of the covariance entries.
    """
    # 1. Pre-processing    
    # Filter SNPs:
    # a) Must not be fixed (mu > 0 and mu < 1)
    # b) Must pass min_maf threshold
    # c) Must not have NaNs (if any exist in row)
    n_snps, n_pops = allele_freq_matrix.shape
    mu = np.nanmean(allele_freq_matrix, axis=1, keepdims=True) # shape (n_snps, 1)

    valid_mask = (
        (mu > min_maf) & 
        (mu < (1.0 - min_maf)) & 
        ~np.isnan(mu)
    ).flatten()

    # Apply filter
    F_valid = allele_freq_matrix[valid_mask]
    mu_valid = mu[valid_mask]
    
    if len(F_valid) == 0:
        raise ValueError("No valid SNPs remaining after filtering.")

    # 2. Per-SNP Normalization (Scaling)
    # Denominator = sqrt( mu * (1-mu) )
    # This standardizes the SNP variance to 1.
    denominator = np.sqrt(mu_valid * (1.0 - mu_valid))
    
    # T is the "Standardized Frequency Matrix" (Z-scores)
    T = (F_valid - mu_valid) / denominator
    
    # 3. Block Splitting
    # Split T into blocks along the SNP axis (axis 0)
    p = n_snps // block_size_snps
    
    if p < 2:
        raise ValueError(f"Not enough SNPs ({n_snps}) to form at least 2 blocks of size {block_size_snps}.")

    limit = p * block_size_snps
    T_truncated = T[:limit, :] # Strictly discard remainder
    blocks = np.split(T_truncated, p)

    # 4. Calculate W for each block (W_ijk)
    W_blocks = []
    
    for block in blocks:
        w_k = (block.T @ block) / block_size_snps
        W_blocks.append(w_k)
        
    W_blocks = np.array(W_blocks) # Shape: (p, n_pops, n_pops)

    # 5. Compute Mean and Standard Error
    W_mean = np.mean(W_blocks, axis=0)
    diff = W_blocks - W_mean 
    sum_sq_diff = np.sum(diff**2, axis=0)
    denominator_se = p * (p - 1)
    
    W_se = np.sqrt(sum_sq_diff / denominator_se)
    
    return W_mean, W_se


def ibd_data_simulated(ts: tskit.TreeSequence, 
            length_bins_cm: list, 
            recombination_rate_cm_per_bp: float):
    """
    Calculates the average IBD length per haploid pair (node pair) 
    for N populations across specified length bins.
    """
    
    pop_nodes = {pop.id: [] for pop in ts.populations()}
    
    for node in ts.nodes():
        # Only look at sample nodes (tips of the tree)
        if node.flags & tskit.NODE_IS_SAMPLE:
            if node.population != tskit.NULL:
                pop_nodes[node.population].append(node.id)

    # Dictionary to hold the sums: stats[(pop_i, pop_j)][bin_index] = total_length
    ibd_sums = {}
    pop_ids = sorted(pop_nodes.keys())
    num_bins = len(length_bins_cm)

    # 2. Iterate through Population Pairs
    # -----------------------------------
    for pop_i, pop_j in combinations_with_replacement(pop_ids, 2):
        ibd_sums[(pop_i, pop_j)] = np.zeros(num_bins)
        
        nodes_i = pop_nodes[pop_i]
        nodes_j = pop_nodes[pop_j]
        
        if not nodes_i or not nodes_j:
            continue

        # 3. Call tskit builtin: ibd_segments
        # -----------------------------------
        if pop_i == pop_j:
            # WITHIN population
            # Returns unique pairs (u, v) where u < v
            ibd_iter = ts.ibd_segments(within=nodes_i, store_pairs=True)
        else:
            # BETWEEN populations
            # Returns pairs (u, v) where u in nodes_i and v in nodes_j
            ibd_iter = ts.ibd_segments(within=nodes_i, between=nodes_j, store_pairs=True)

        # 4. Process Segments
        # -------------------
        for segment in ibd_iter:
            # Convert bp to cM
            length_cm = (segment.right - segment.left) * recombination_rate_cm_per_bp
            
            # Find the correct bin
            for bin_idx, (min_cm, max_cm) in enumerate(length_bins_cm):
                if min_cm <= length_cm < max_cm:
                    # Accumulate the LENGTH (as requested)
                    ibd_sums[(pop_i, pop_j)][bin_idx] += length_cm
                    break

    # 5. Calculate Averages (Normalization)
    # -------------------------------------
    results = []
    
    for (pop_i, pop_j), sums in ibd_sums.items():
        # Get Haploid Counts (H)
        H_i = len(pop_nodes[pop_i])
        H_j = len(pop_nodes[pop_j])
        
        # Calculate denominator based on 40 * 40 logic
        total_pairs = H_i * H_j
        
        # ---------------------------------------------------------
        # CRITICAL SYMMETRY STEP
        # ---------------------------------------------------------
        # Case 1: Within Population (pop_i == pop_j)
        # ibd_segments returns unique pairs (Triangle).
        # To average over the full Square Matrix (Ni * Ni), we must
        # double the sum of the off-diagonal elements.
        # (We assume self-IBD is excluded or handled as 0 in this specific statistic)
        if pop_i == pop_j:
            final_sum = sums * 2
        else:
            # Case 2: Between Populations
            # ibd_segments returns all links between Group A and Group B.
            # No doubling needed.
            final_sum = sums

        if total_pairs > 0:
            averages = final_sum / total_pairs
        else:
            averages = [0] * num_bins
            
        for bin_idx, avg in enumerate(averages):
            results.append({
                "pop_1": pop_i,
                "pop_2": pop_j,
                "bin_min_cm": length_bins_cm[bin_idx][0],
                "bin_max_cm": length_bins_cm[bin_idx][1],
                "avg_ibd_length": avg,
                "haploids_1": H_i,
                "haploids_2": H_j
            })
            
    return results