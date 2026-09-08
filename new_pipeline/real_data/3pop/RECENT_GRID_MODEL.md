# Recent Ne grid model

The `Ngrid` variants add two fixed, absolute recent epochs to every sampled
population:

| epoch | generations |
|---|---:|
| 1 | 0-5 |
| 2 | 5-10 |

After generation 10, each leaf uses its ordinary branch-level Ne until the first
demographic event. Internal and ancestry-source branches retain the existing
tree-structured Nsmooth model.

## Event-time convention

The first sampled time parameter is the excess beyond generation 10:

```text
first event time = 10 + first excess time
```

The excess has the existing `exponential(0.01)` prior and a one-generation lower
bound, so the first event is at least generation 11. This is the exponential
prior conditional on the first event occurring after the fixed grid. Later time
parameters remain gaps between consecutive events.

This permits inferred events throughout the 10-25 generation window without
allowing an event to cross a fixed grid boundary. Supporting events inside
0-10 would require event-dependent lineage grids and a discontinuous change in
which Ne parameters are active when an event crosses generation 5 or 10.

## Smoothness prior

For a five-generation step, the existing random-walk scale is reused:

```text
log Ne[5-10] = log Ne[post-10 leaf branch]
               + tau * sqrt(5 / 100) * z[5-10]

log Ne[0-5]  = log Ne[5-10]
               + tau * sqrt(5 / 100) * z[0-5]

z ~ standard Normal
```

Thus no new smoothness hyperparameter is introduced. Shared-Ne models use one
recent trajectory for SNP and IBD; separate-Ne models use independent recent
trajectories and their existing `tau_ibd` and `tau_snp` parameters.

## Models

- `mixed_model_Ngrid_poisson.stan`
- `mixed_model_Ngrid_normal.stan`
- `mixed_model_Ngrid_poisson_separate_ne.stan`
- `mixed_model_Ngrid_normal_separate_ne.stan`

Run all four grid variants over all 21 topologies with:

```bash
/opt/anaconda3/envs/stan_env/bin/python 3pop/fit_four_models.py \
  --tag eas_ibs_tsi --pops EAS IBS TSI \
  --variants poisson_grid_shared_ne normal_grid_shared_ne \
             poisson_grid_separate_ne normal_grid_separate_ne \
  --map-starts 12 --pathfinder-modes 3 --map-iter 1200 \
  --draws 4000 --paths 8
```

All 21 topologies in `eas_ibs_tsi` were fitted for all four grid variants with
12 dispersed MAP starts, three promoted modes, 4,000 draws, and eight
Pathfinder paths. The combined outputs are
in `eas_ibs_tsi/comparison/report.md` and
`eas_ibs_tsi/comparison/eight_model_topology_table.csv`.
