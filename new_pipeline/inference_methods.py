import numpy as np
import tskit
from collections import defaultdict
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


def prepare_snp_blocks(mts, dem, block_size_cm, cm_per_unit, cutoff_time=None, min_maf=0.05):
    """
    Extract per-SNP normalized frequency vectors, grouped by genomic block.

    Each genomic block (50 cM, matching IBD) stores a list of per-SNP
    t-vectors.  At resampling time the chosen blocks' SNPs are gathered
    and TreeMix SE is computed from 500-SNP sub-blocks.

    Parameters
    ----------
    mts : tskit.TreeSequence
    dem : DemographicTopology
    block_size_cm : float
        Genomic block size in cM (must match IBD block size).
    cm_per_unit : float
    cutoff_time : float, optional
    min_maf : float

    Returns
    -------
    snp_block_t : list of np.ndarray
        snp_block_t[k] has shape (n_snps_in_block_k, n_pops).
    n_blocks : int
    """
    leaves = dem.initial_leaves
    n_pops = len(leaves)

    pop_sample_ids = {}
    for pop_name in leaves:
        for p in mts.populations():
            if p.metadata.get('name') == pop_name:
                pop_sample_ids[pop_name] = mts.samples(population=p.id)
                break

    genome_length_bp = mts.sequence_length
    genome_length_cm = genome_length_bp * cm_per_unit
    block_size_bp = block_size_cm / cm_per_unit
    n_blocks = int(np.floor(genome_length_cm / block_size_cm))
    used_length_bp = n_blocks * block_size_bp

    if n_blocks == 0:
        raise ValueError(
            f"Genome ({genome_length_cm:.1f} cM) shorter than one block ({block_size_cm} cM)."
        )

    block_dev_lists = [[] for _ in range(n_blocks)]
    block_het_lists = [[] for _ in range(n_blocks)]
    n_filtered = 0
    n_kept = 0

    for variant in mts.variants():
        if cutoff_time is not None:
            mut_node = variant.site.mutations[0].node
            if mts.node(mut_node).time < cutoff_time:
                n_filtered += 1
                continue

        pos = variant.site.position
        if pos >= used_length_bp:
            continue
        block_idx = int(pos // block_size_bp)
        if block_idx >= n_blocks:
            continue

        genotypes = variant.genotypes
        freqs = np.array([
            np.mean(genotypes[pop_sample_ids[pop]])
            for pop in leaves
        ])

        mu = np.mean(freqs)
        if mu <= min_maf or mu >= (1.0 - min_maf):
            n_filtered += 1
            continue

        # Store unnormalized deviation and heterozygosity separately
        # (TreeMix pooled normalization: ratio of sums, not sum of ratios)
        dev = freqs - mu
        het = mu * (1.0 - mu)
        block_dev_lists[block_idx].append(dev)
        block_het_lists[block_idx].append(het)
        n_kept += 1

    snp_blocks = []
    for devs, hets in zip(block_dev_lists, block_het_lists):
        if len(devs) > 0:
            snp_blocks.append((np.array(devs), np.array(hets)))
        else:
            snp_blocks.append((np.empty((0, n_pops)), np.empty(0)))

    print(f"SNP blocks: {n_kept} SNPs kept across {n_blocks} blocks "
          f"({n_filtered} filtered), avg {n_kept / max(n_blocks, 1):.0f} SNPs/block")

    return snp_blocks, n_blocks


def pool_snp_blocks(snp_block_list):
    """
    Concatenate per-block SNP data from multiple simulations.

    Parameters
    ----------
    snp_block_list : list of list[tuple(np.ndarray, np.ndarray)]
        Each element is a simulation's snp_blocks (list of (dev, het) tuples).

    Returns
    -------
    pooled : list of tuple(np.ndarray, np.ndarray)
        Flat list of per-genomic-block (dev, het) tuples.
    """
    pooled = []
    for block_list in snp_block_list:
        pooled.extend(block_list)
    total_snps = sum(blk[0].shape[0] for blk in pooled)
    print(f"Pooled SNP blocks: {len(pooled)} genomic blocks, "
          f"{total_snps} total SNPs")
    return pooled


def resample_snp_covariance(snp_blocks, chosen_blocks,
                            n_haploid_per_pop=None,
                            se_block_size=500, se_floor=1e-8):
    """
    Compute w_hat and w_se from resampled genomic blocks (TreeMix pooled approach).

    Uses the pooled normalization matching the TreeMix software:
        w_ij = sum_a (x_i - x_bar)(x_j - x_bar) / sum_a x_bar(1 - x_bar)
    This is a ratio of sums (SNPs weighted by heterozygosity),
    NOT a sum of ratios (equal weight per SNP).

    Parameters
    ----------
    snp_blocks : list of tuple(np.ndarray, np.ndarray)
        Per-genomic-block (dev, het) tuples from pool_snp_blocks.
        dev[k] has shape (n_snps_k, n_pops): frequency deviations (x - x_bar).
        het[k] has shape (n_snps_k,): heterozygosities x_bar*(1-x_bar).
    chosen_blocks : np.ndarray of int
        Genomic block indices selected (same as IBD resampling).
    n_haploid_per_pop : array-like (n_pops,), optional
        For TreeMix finite-sample correction.
    se_block_size : int
        Number of SNPs per SE block (TreeMix window_size, default 500).
    se_floor : float
        Minimum SE.

    Returns
    -------
    w_hat : np.ndarray (n_pops, n_pops)
    w_se  : np.ndarray (n_pops, n_pops)
    """
    # 1. Gather deviations and heterozygosities from chosen genomic blocks
    dev_arrays = [snp_blocks[k][0] for k in chosen_blocks]
    het_arrays = [snp_blocks[k][1] for k in chosen_blocks]
    D = np.concatenate(dev_arrays, axis=0)   # (J, n_pops)
    H = np.concatenate(het_arrays, axis=0)   # (J,)
    J, n_pops = D.shape

    if J == 0:
        return np.zeros((n_pops, n_pops)), np.full((n_pops, n_pops), se_floor)

    # 2. Pooled covariance: w_hat = sum(d_a d_a') / sum(h_a)
    #    (already double-centered because d_a = freqs - mean(freqs)
    #     is zero-sum across populations by construction)
    numer = D.T @ D          # (n_pops, n_pops)
    denom = np.sum(H)        # scalar
    w_hat = numer / denom

    # 3. SE from sub-blocks (each block uses its own pooled covariance)
    B = J // se_block_size
    if B >= 2:
        D_trunc = D[:B * se_block_size]
        H_trunc = H[:B * se_block_size]
        D_chunks = D_trunc.reshape(B, se_block_size, n_pops)
        H_chunks = H_trunc.reshape(B, se_block_size)

        numer_blocks = np.einsum('bsp,bsq->bpq', D_chunks, D_chunks)  # (B, n_pops, n_pops)
        denom_blocks = H_chunks.sum(axis=1)                            # (B,)
        w_blocks = numer_blocks / denom_blocks[:, np.newaxis, np.newaxis]

        w_mean = w_blocks.mean(axis=0)
        diff = w_blocks - w_mean[np.newaxis]
        sum_sq = np.sum(diff ** 2, axis=0)
        w_se = np.sqrt(np.maximum(sum_sq / (B * (B - 1)), se_floor ** 2))
    else:
        w_se = np.full((n_pops, n_pops), se_floor)

    # 4. TreeMix finite-sample correction (Text S1 Eqs. 7-8):
    #    With pooled normalization, diagonal bias is still 1/n_i
    #    (sampling variance sum_a p(1-p)/n_i cancels with denom sum_a p(1-p)).
    #    Double-center the correction to stay in centered space.
    if n_haploid_per_pop is not None:
        bias = 1.0 / np.asarray(n_haploid_per_pop, dtype=float)
        noise_raw = np.diag(bias)
        row_means = noise_raw.mean(axis=1, keepdims=True)
        grand_mean = noise_raw.mean()
        noise_centered = noise_raw - row_means - row_means.T + grand_mean
        w_hat = w_hat - noise_centered

    return w_hat, w_se


