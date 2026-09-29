"""Stan data dictionaries for the IBD-only, SNP-only and mixed models.

All three share a topology block (events as migration matrices, admixture map,
index arrays) and the leaf-pair list; each adds its own observations.  The
Nsmooth models also need the tree the log-Ne random walk runs on (`smooth_data`).
"""
from __future__ import annotations

import numpy as np


def topology_data(dem):
    """Topology block shared by every model."""
    matrices, _, _, admix_id = dem.get_topology_matrix_representation()
    n_events = len(dem.ordered_events)
    admixture_map = [[k, admix_id[k]["child"]] + admix_id[k]["parents"] for k in admix_id]
    admixture_indices = [k + 2 for k in admix_id]
    fixed_shifted = [k for k in range(2, n_events + 2) if k not in admixture_indices]
    L = len(dem.initial_leaves)
    pair_i, pair_j = (list(x) for x in zip(*[(i, j) for i in range(1, L + 1) for j in range(i, L + 1)]))
    return dict(n_leaves=L, n_nodes=len(dem.nodes), n_events=n_events, n_admixture=dem.n_admix,
                migration_matrices=matrices,
                admixture_map=np.asarray(admixture_map, dtype=int).reshape(-1, 4),
                admixture_indices=admixture_indices, fixed_indices_shifted=fixed_shifted,
                fixed_indices=[k - 1 for k in fixed_shifted],
                n_leaf_pairs=len(pair_i), pair_i=pair_i, pair_j=pair_j)


def smooth_data(dem):
    """Tree for the Nsmooth log-Ne random walk (1-based; parents have higher index).

    ne_parent       node each branch flows into going back in time (0 = root/admix source)
    ne_admix_idx    admixture index for a branch created as an admixture source (else 0)
    ne_start_event  event that creates the branch (0 = leaf, starts at t = 0)
    """
    idx = {name: i for i, name in enumerate(dem.nodes)}
    parent, admix, start = ([0] * len(idx) for _ in range(3))
    a = 0
    for k, ev in enumerate(dem.ordered_events, 1):
        if ev["type"] == "MERGE":
            start[idx[ev["parent"]]] = k
            for child in ev["children"]:
                parent[idx[child]] = idx[ev["parent"]] + 1
        else:
            a += 1
            admix[idx[ev["child"]]] = a
            for name in ev["parents"]:
                start[idx[name]] = k
    return dict(ne_parent=parent, ne_admix_idx=admix, ne_start_event=start)


def ibd_data(edges, obs, n_haploid, T_max, se_floor=1e-8):
    """IBD observations.  obs is `blocks.aggregate` output (ibd_hat, ibd_se, cm, ibd_count);
    `stan_data` keeps only what the model declares (counts for Poisson, hat/se for normal)."""
    edges = np.asarray(edges, float)
    out = dict(n_bins=len(edges) - 1, bin_length=[[lo, hi] for lo, hi in zip(edges[:-1], edges[1:])],
               T_max=T_max, cm=obs["cm"], n_samples=[int(n) for n in n_haploid],
               ibd_hat=np.asarray(obs["ibd_hat"]).tolist(),
               ibd_se=np.maximum(np.asarray(obs["ibd_se"]), se_floor).tolist(),
               ibd_count=np.asarray(obs["ibd_count"], dtype=int).tolist())
    return out


def snp_data(obs):
    return dict(w_hat=np.asarray(obs["w_hat"]).tolist(), w_se=np.asarray(obs["w_se"]).tolist())


def stan_data(dem, obs, spec, edges=None, n_haploid=None, T_max=1e5):
    """Complete data for model `spec` (a `models.ModelSpec`) fitted on topology `dem`:
    exactly the variables its `data {}` block declares, no more."""
    data = topology_data(dem)
    if spec.data in ("ibd", "mixed"):
        data.update(ibd_data(edges, obs, n_haploid, T_max))
    if spec.data in ("snp", "mixed"):
        data.update(snp_data(obs))
    if spec.ne == "smooth":
        data.update(smooth_data(dem))
    names = spec.data_names
    missing = names - set(data)
    if missing:
        raise ValueError(f"{spec.file} needs data this builder does not provide: {sorted(missing)}")
    return {k: v for k, v in data.items() if k in names}
