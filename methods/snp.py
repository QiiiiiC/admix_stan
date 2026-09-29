"""SNP allele-frequency covariance (TreeMix-style), with a block jackknife SE.

Per genome block we keep additive sums over SNP sub-blocks, so any subsample of
blocks can be aggregated without touching genotypes again:

  numer[k] = sum_snps outer(f - fbar, f - fbar)     (L x L)
  denom[k] = sum_snps fbar (1 - fbar)

W = sum numer / sum denom (ratio of sums), minus the finite-sample correction
C diag(1/n_i) C, where C is the L x L centering matrix.  The SE is a
delete-one-sub-block jackknife over every sub-block of the chosen blocks.
"""
from __future__ import annotations

import numpy as np


def snp_summaries(ts, pops, block_cm, sub_block_cm, recombination_rate, min_maf=0.05,
                  min_mutation_time=None):
    """Per-sub-block numer/denom/site counts for one simulated block.

    Biallelic sites only.  min_mutation_time: keep only sites carrying a single
    mutation on a node at least this old (ascertain SNPs that predate the
    graph's root, as in an outgroup-ascertained panel); None keeps every site.
    The MAF filter is on the mean frequency across populations.
    """
    pop_id = {p.metadata["name"]: p.id for p in ts.populations()}
    column = {int(u): j for j, u in enumerate(ts.samples())}
    indices = [[column[int(u)] for u in ts.samples(population=pop_id[p])] for p in pops]
    L = len(pops)
    n = round(block_cm / sub_block_cm)
    numer = np.zeros((n, L, L)); denom = np.zeros(n); sites = np.zeros(n, np.int64)
    for var in ts.variants():
        if len(var.alleles) != 2 or np.any(var.genotypes < 0):
            continue
        if min_mutation_time is not None and (len(var.site.mutations) != 1 or
                ts.node(var.site.mutations[0].node).time < min_mutation_time):
            continue
        freq = np.asarray([np.mean(var.genotypes[ix]) for ix in indices])
        mean = freq.mean()
        if not min_maf < mean < 1 - min_maf:
            continue
        k = min(n - 1, int(var.site.position * 100 * recombination_rate / sub_block_cm))
        dev = freq - mean
        numer[k] += np.outer(dev, dev)
        denom[k] += mean * (1 - mean)
        sites[k] += 1
    return dict(snp_numer=numer, snp_denom=denom, snp_sites=sites)


def snp_covariance(numer, denom, n_haploid):
    """(W, SE) from stacked sub-block sums.  n_haploid: haploid sample size per population."""
    n_haploid = np.asarray(n_haploid, float)
    L = len(n_haploid)
    numer = np.asarray(numer).reshape(-1, L, L)
    denom = np.asarray(denom).reshape(-1)
    if len(denom) < 2 or np.sum(denom) <= 0 or np.any(denom.sum() - denom <= 0):
        raise ValueError("Not enough SNP-bearing blocks to estimate uncertainty")
    center = np.eye(L) - np.ones((L, L)) / L
    correction = center @ np.diag(1 / n_haploid) @ center
    w = numer.sum(axis=0) / denom.sum() - correction
    leave = (numer.sum(axis=0) - numer) / (denom.sum() - denom)[:, None, None] - correction
    se = np.sqrt((len(denom) - 1) / len(denom) * np.sum((leave - leave.mean(axis=0))**2, axis=0))
    return w, np.maximum(se, 1e-8)
