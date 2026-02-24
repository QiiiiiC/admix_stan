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


