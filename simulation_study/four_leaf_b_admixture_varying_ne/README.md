# Four leaves, one admixture, a different Ne on every branch

Parameter recovery only.  The generating graph is fitted as the one and only
candidate, and the question is which of its 16 parameters each model recovers,
and from what genome length.  The parameters are 5 event times, 1 admixture
fraction and 10 branch sizes.  This follows `admix_thru_b_plus_ancient`, where
one Ne was shared by every branch (the Nfixed models).

## Generating history (`summary/figures/topology_true.png`)

Populations `a, b, c, d`, 20 diploid samples each.  Backwards in time:

| Generations | Event | Parameter |
|---:|---|---|
| 20 | b splits into b1 (0.7) and b2 (0.3) | `t_b_admix`, `f_b1` |
| 60 | a and b1 merge → ab | `t_ab` |
| 100 | b2 and c merge → cb | `t_cb` |
| 200 | cb and d merge → cbd | `t_cbd` |
| 400 | ab and cbd merge → root | `t_root` |

Haploid Ne per branch (`Ne_<branch>`), a 4-fold range with neighbours differing
by up to 4×:

| a | b | c | d | b1 | b2 | ab | cb | cbd | root |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20,000 | 5,000 | 12,000 | 8,000 | 16,000 | 6,000 | 12,000 | 10,000 | 15,000 | 15,000 |

b is small and admixed from a large source (b1) and a small one (b2).  The
events span 20–400 generations, so IBD (segments ≥ 1 cM) has leverage on most
of the graph.  It has almost none on the root branch: coalescence 400+
generations ago leaves segments of ~0.1 cM, so **`Ne_root` is expected to be
prior-driven**.

## Models (all with the Poisson IBD likelihood on segment counts)

| Key | Model | Ne |
|---|---|---|
| `ibd` | `ibd_Nvarying_poisson` | an independent Ne per branch, log-normal around a shared level (`mu_log`, `sigma_log`) |
| `mixed` | `mixed_Nvarying_poisson` | the same, with the SNP covariance term |
| `mixed_smooth` | `mixed_Nsmooth_poisson` | Ne as a random walk along the tree (the latest model), for comparison |

There is **no SNP-only model**.  SNP drift on a branch is duration / Ne, so SNP
alone identifies only the ratio.  This was checked in
`admix_thru_b_plus_ancient`: SNP-only Ne and root time are both 9–13% low, and
their ratio is right to < 1%.

Nsmooth's prior centres each branch on its parent with step sd
`tau · sqrt(dt / 100)` (`tau ~ half-normal(0.3)`).  A 4× jump on a short branch
is therefore a many-sigma excursion, and fitting it measures how much that
smoothing shrinks a truth like this one.  The `*_Nvarying_poisson` models are
derived from the normal-likelihood Nvarying files by
`methods/stan/make_poisson.py`.

## Design

Lengths 50–1000 cM (the same grid as before), 100 subsamples each, drawn
from one pool of 100 × 25 cM blocks.  Both true IBD and hap-IBD are used, with
the same bins and SNP ascertainment as `admix_thru_b_plus_ancient`.  There are
three dispersed, truth-free Pathfinder starts per fit.  That makes
8 × 100 × 3 × 2 = **4,800 fits**.

## Outputs (`summary/`)

- `report.md`: per IBD source, at the shortest and longest length, the median
  relative error and the 95%-interval coverage of every parameter, by model.
- `recovery.csv`: the same for every length.  `parameter_estimates.csv` has
  every fit × parameter; `fits.csv` has the diagnostics.
- `figures/parameters_*`: relative error vs length, 16 panels.
  `figures/coverage_*`: coverage vs length.  `figures/topology_true`: the
  graph with every branch named and its Ne.

## Running

```sh
python execute_full.py --jobs 8
python -B -m unittest test_design -v
```

Stages and resumption are as in `admix_thru_b_plus_ancient`.  `calibrate`
runs NUTS on chosen lengths (`--calibration-cm`) as a check on Pathfinder's
intervals.

## Caveats

Intervals come from Pathfinder importance weights, whose ESS is often low, so
coverage is indicative, not calibrated.  Subsamples overlap.  The Nvarying prior
shares one level across branches, so it also shrinks, only less than Nsmooth
does.

## Rerun through the growth-design pipeline (`runs/pipeline`, `summary/pipeline`)

The same design (same trees, same subsamples) refitted with the pipeline in
`../four_leaf_b_admixture_growth_ne`, so the two designs are directly comparable:

- `mixed_Nsmooth_poisson` and `ibd_Nsmooth_poisson`, on T_true and the no-admixture
  tree, ranked by best-path ELBO;
- 22 bins from 1.0 cM, with IBD re-extracted from the saved tree sequences;
- NUTS + bridge sampling on 10 subsamples per length for the mixed model;
- spectrum fit at 1000 cM and ELBO vs log Z.

Config: `pipeline_config.json`.  See that folder's README for the outputs.
