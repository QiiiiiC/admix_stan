# Four leaves, one admixture, a different Ne on every branch

The companion of `../four_leaf_b_admixture_growth_ne`: the same graph (b admixed at
20 gen, 0.7 from the a side, merges at 60 / 100 / 200 / 400), bins, models, NUTS
comparison and figures, run by that folder's code.  Only the branch sizes differ:
an explicit table, smaller and unordered in time, instead of growth.

| a | b | c | d | b1 | b2 | ab | cb | cbd | root |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20,000 | 5,000 | 12,000 | 8,000 | 16,000 | 6,000 | 12,000 | 10,000 | 15,000 | 15,000 |

(haploid; `config.json`, seed 20260929, 100 × 25 cM blocks, 20 diploids per leaf.)

## What is fitted

| | |
|---|---|
| Models | `mixed` = `mixed_Nsmooth_poisson`, `ibd` = `ibd_Nsmooth_poisson`.  Both Nsmooth, Poisson IBD likelihood. |
| Graphs | `T_true`; `T_null` = ((a,b),(c,d)), b wholly on its major source side. |
| IBD sources | `true` and `hapibd` (hap-IBD on the simulated genotypes, min-seed 1.0, min-output 1.0 cM). |
| Pathfinder | both sources × both models × both graphs × 100 subsamples × 8 lengths (50–1000 cM) = 6,400 fits; graph score = best-path ELBO. |
| NUTS + bridge | true IBD, mixed × both graphs × 10 subsamples × 8 lengths = 160 fits. |
| Bins | 22, from 1.0 cM: 0.5 cM steps to 10 cM, then 10–12, 12–15, 15–20, 20–25. |

## Outputs (`summary/`)

Same as the growth design: `topology_rank`, `parameters_<source>`,
`spectrum_1000cM_<source>`, `residuals_1000cM_<source>`, `hapibd_detection`,
`elbo_vs_logz`, `topology_true`; tables in `report.md` and the CSVs.

## Run

From `../four_leaf_b_admixture_growth_ne`:

```sh
python execute_full.py --config ../four_leaf_b_admixture_varying_ne/config.json \
    --out ../four_leaf_b_admixture_varying_ne/runs/default \
    --publish ../four_leaf_b_admixture_varying_ne/summary --stages simulate fit nuts summarize
```

`runs/default/MIGRATION.json` records how the run reached this layout (it was first
fitted from tree sequences simulated by this design's original runner, with the same
seeds; re-simulating a block with the pipeline reproduces them bit for bit).
