"""Borrowed data: re-extract true IBD (in this folder's bins) and SNP sums from the
tree sequences of another study, and aggregate subsamples of whole blocks."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from methods import blocks, ibd, snp
from methods.utils import save_json


def load_config(path=HERE / "config.json"):
    return json.loads(Path(path).read_text())


def n_haploid(c):
    return [2 * c["diploid_samples"]] * len(c["leaves"])


def extract(c, out, block_slice=None):
    import tskit
    folder = out / "pool_000"; folder.mkdir(parents=True, exist_ok=True)
    source = ROOT / c["trees_from"]
    save_json(folder / "provenance.json", {"trees_from": str(source), "tskit": tskit.__version__})
    for k in range(c["blocks"]):
        if block_slice and k % block_slice[1] != block_slice[0]:
            continue
        cache = folder / f"block_{k:03d}.npz"
        if cache.exists() and snp.is_current(np.load(cache)):
            continue
        ts = tskit.load(source / f"block_{k:03d}.trees")
        values = snp.snp_summaries(ts, c["leaves"], c["block_cm"], c["snp_se_block_cm"], c["recombination_rate"],
                                   c["snp_min_maf"], c["snp_min_mutation_time"])
        values["true_count"], values["true_length"] = ibd.true_ibd(ts, c["leaves"], np.asarray(c["bin_edges_cm"]),
                                                                   c["recombination_rate"])
        tmp = cache.with_suffix(".tmp.npz"); np.savez_compressed(tmp, **values); tmp.replace(cache)
        print(f"block {k}: {int(values['true_count'].sum())} segments", flush=True)


def subsample(c, out, cm, rep):
    pool = blocks.load_blocks(out / "pool_000", c["blocks"], ["true"])
    n = round(cm / c["block_cm"])
    chosen = blocks.selections(c["seed"], 0, c["blocks"], n, c["replicates"], int(cm))[rep]
    return blocks.aggregate(pool, chosen, "true", c["block_cm"], n_haploid(c)), chosen


def stan_data(c, sk, obs):
    """Everything the models in stan/ declare, except model-specific entries."""
    e = np.asarray(c["bin_edges_cm"], float)
    return dict(L=sk.L, n_merge=sk.n_merge, n_pairs=len(sk.pairs), pairs=[[a + 1, b + 1] for a, b in sk.pairs],
                n_bins=len(e) - 1, bin_length=[[lo, hi] for lo, hi in zip(e[:-1], e[1:])], T_max=c["T_max"],
                cm=float(obs["cm"]), n_samples=n_haploid(c),
                ibd_count=np.asarray(obs["ibd_count"], int).tolist(),
                w_hat=np.asarray(obs["w_hat"]).tolist(), w_se=np.asarray(obs["w_se"]).tolist(), use_snp=c["use_snp"])
