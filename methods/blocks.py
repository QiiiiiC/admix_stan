"""Block pool -> observed data for one subsample.

A study simulates a pool of independent genome blocks once, caches each
block's IBD counts/lengths and SNP sums (`block_XXX.npz`), and then forms every
observation by choosing blocks without replacement.  Varying the number of
blocks chosen is how genome length is varied.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .ibd import pair_counts
from .snp import snp_covariance


def selections(seed, pool, n_blocks, n_chosen, replicates, *key):
    """`replicates` sorted subsets of n_chosen out of n_blocks block ids.

    Seeded by (seed, pool, 987, *key): the same ids serve every IBD source,
    model and scenario.  Pass the genome length (or anything else that should
    draw a fresh set) as extra key entries; with no key the stream is the one
    growth_asymmetric_loop used.
    """
    rng = np.random.default_rng(np.random.SeedSequence([seed, pool, 987, *key]))
    return [sorted(rng.choice(n_blocks, n_chosen, replace=False).tolist()) for _ in range(replicates)]


def load_blocks(folder, n_blocks, sources):
    blocks = []
    for k in range(n_blocks):
        path = Path(folder) / f"block_{k:03d}.npz"
        with np.load(path) as z:
            block = dict(z)
        if any(s + "_count" not in block for s in sources):
            raise ValueError(f"Missing IBD source in {path}; finish simulate first")
        blocks.append(block)
    return blocks


def aggregate(blocks, chosen, source, block_cm, n_haploid):
    """Observed data for one subsample and one IBD source.

    ibd_count  summed segment counts per bin and pair (what Poisson models use)
    ibd_hat    IBD fraction: summed length / (genome cM * pair exposure)
    ibd_se     SE of ibd_hat across the chosen blocks (normal-likelihood models)
    w_hat/w_se SNP covariance and its sub-block jackknife SE
    cm         genome length of the subsample
    """
    exposure = pair_counts(n_haploid)
    count = np.sum([blocks[i][source + "_count"] for i in chosen], axis=0)
    lengths = np.asarray([blocks[i][source + "_length"] for i in chosen])
    cm = len(chosen) * block_cm
    fraction = lengths.sum(axis=0) / (cm * exposure)
    per_block = lengths / (block_cm * exposure)
    se = np.maximum(per_block.std(axis=0, ddof=1) / np.sqrt(len(chosen)), 1e-8)
    w, ws = snp_covariance([blocks[i]["snp_numer"] for i in chosen],
                           [blocks[i]["snp_denom"] for i in chosen], n_haploid)
    return dict(ibd_count=count, ibd_hat=fraction, ibd_se=se, w_hat=w, w_se=ws, cm=cm)
