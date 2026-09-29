# Four leaves, one admixture, growth-like Ne: topology, recovery, ELBO vs log Z

Same graph as `four_leaf_b_admixture_varying_ne` (b admixed at 20 gen, 0.7 from the
a side, merges at 60 / 100 / 200 / 400), but the recent branches are **large and
different**.  The sizes mimic exponential growth while staying constant along each
branch:

- every lineage grows log-linearly from 10,000 at the root (400 gen) to its own
  present-day size (a 150k, b 50k, c 100k, d 30k);
- a branch shared by several leaves follows the geometric mean of their present sizes;
- each branch takes the curve's value at its time midpoint (`study.branch_ne`).

| a | b | c | d | b1 | b2 | ab | cb | cbd | root |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 122,429 | 48,028 | 74,989 | 22,795 | 42,567 | 39,276 | 25,029 | 33,957 | 15,182 | 10,000 |

## What is fitted

| | |
|---|---|
| Models | `mixed` = `mixed_Nsmooth_poisson`, `ibd` = `ibd_Nsmooth_poisson` (IBD only).  Both Nsmooth, Poisson IBD likelihood. |
| Graphs | `T_true`; `T_null` = ((a,b),(c,d)), b wholly on its major source side.  A graph is an event order. |
| IBD sources | `true` (segments read off the genealogies) and `hapibd` (hap-IBD on the simulated genotypes, min-seed 1.0, min-output 1.0 cM). |
| Pathfinder | both sources × both models × both graphs × 100 subsamples × 8 lengths (50–1000 cM) = 6,400 fits.  Graph score = **best-path ELBO**. |
| NUTS + bridge | true IBD, mixed × both graphs × 10 subsamples × 8 lengths = 160 fits: exact log Z (`methods.evidence`), parameters, posterior predictive. |
| Bins | 22, from 1.0 cM: 0.5 cM steps to 10 cM, then 10–12, 12–15, 15–20, 20–25. |

## Outputs (`summary/`)

- `topology_rank`: ΔELBO(T_true − T_null) per model vs length, and the identification
  rate (95% Wilson interval); one row per IBD source.
- `parameters_true` / `parameters_hapibd`: relative error of all 16 T_true
  parameters, from Pathfinder (both models, 100 subsamples) and, for true IBD, NUTS
  (mixed, 10 subsamples).
- `spectrum_1000cM_<source>` / `residuals_1000cM_<source>`: observed counts in every
  bin and pair across all 100 subsamples at 1000 cM, with both models' fitted
  expectations, and the Pearson residuals.
- `hapibd_detection`: hap-IBD segment counts over true counts per bin, all blocks.
- `elbo_vs_logz`: for the mixed model, log Z − best-path ELBO per graph, paired
  ΔlogZ vs ΔELBO, and identification by ELBO vs by log Z.
- `topology_true`: the generating graph with every branch's Ne.
- Tables: `report.md`, `fits.csv`, `topology.csv`, `parameter_estimates.csv`,
  `hapibd_detection.csv`, `elbo_vs_logz.csv`.

Results live in `runs/<out>/pool_000/<stage>/cm_<L>/rep_<r>/<source>/<model>/<graph>/`.
A block cache is complete once it holds every configured source; `simulate` adds only
the missing ones from the saved tree sequences.

## One pipeline, two designs

The code here also reruns the varying-Ne design
(`../four_leaf_b_admixture_varying_ne/pipeline_config.json`).  Its config names
`trees_from`, so `simulate` re-reads that study's saved tree sequences and
re-extracts IBD in these 22 bins, rather than simulating.  Its subsamples are
identical to the original study's.

```sh
python execute_full.py --stages simulate fit summarize     # Pathfinder part
python execute_full.py --stages nuts summarize              # NUTS part
python -B -m unittest test_design -v
```
