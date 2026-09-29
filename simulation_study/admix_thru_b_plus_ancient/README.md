# Recent admix-through-b loop + ancient admixture (Nfixed)

The case in which the **mixed model beats both single-data models**: a recent
event only IBD can see and an ancient event only SNP covariance can see, in one
graph.  This reruns `new_pipeline/Nfixed/topology_admix_thru_b_plus_ancient`
(the APM report's two-admixture experiment, `sim2`) on the shared `methods/`
pipeline, across genome lengths, with true and hap-IBD-called segments.

## Generating history

Populations `a, b, c, d`, 20 diploid (40 haploid) samples each, haploid Ne
10,000 on every branch.  Backwards in time:

| Generations | Event |
|---:|---|
| 20 | a splits into aP1 (0.5) and aP2 (0.5) — recent admixture |
| 60 | aP2 merges into b (bP) — the donor branch runs through b |
| 100 | bP and aP1 merge (ab) — loop closes |
| 500 | c splits into cP1 (0.7) and cP2 (0.3) — ancient admixture |
| 700 | ab and cP1 merge (left) |
| 900 | d and cP2 merge (right) |
| 1000 | root |

## What is fitted

Three candidates, each in its generating event order: `T_true`; `T_alt1`,
without the recent loop (a and b merge directly); `T_alt2`, without the ancient
admixture (c joins (ab, d) as a tree).  Each is fitted with
`ibd_Nfixed_poisson`, `snp_Nfixed` and `mixed_Nfixed_poisson`, using the
normalised copies from `methods.models`.  All three have one shared
`effective_N` with a gamma prior.  So there is a single Ne parameter for every branch, which is
the right model class here (the truth is 10,000 on every branch), and the `Ne`
panel of `parameters_*` is that one parameter, not an average over branches.
No model is given the true Ne; every start is drawn without it.  IBD enters through the Poisson likelihood
on raw segment counts, the same one the latest studies use, so no IBD
standard error is involved.

The two contrasts are
**ΔELBO(T_true − T_alt1)**, which asks whether the recent loop is kept, and
**ΔELBO(T_true − T_alt2)**, which asks whether the ancient admixture is kept.
The prediction is that IBD-only wins only the first, SNP-only only the second,
and Mixed both.

Genome lengths 50, 100, 150, 200, 300, 500, 750, 1000 cM (the original grid).
One pool of 100 × 25 cM blocks (2,500 cM) is simulated once; each length
draws 100 subsamples of blocks without replacement, with its own seed key.
Blocks are 25 cM (the original used 50).  Bins are the original ones from 1 cM, with
the open top bin closed at the block length (20–25 cM).  SNPs are ascertained
as mutations older than the root, MAF ≥ 0.05, with 5 cM jackknife sub-blocks.

Every fit uses three dispersed, truth-independent Pathfinder starts (4 paths ×
1,000 draws each), and ranking is by ELBO.  logZ, ESS and Pareto k are
recorded per fit.  SNP-only never reads IBD, so it is fitted once per
subsample (`source = none`) and appears in both the true-IBD and hap-IBD
summaries.  That makes 8 × 100 × 3 × (2 + 2 + 1) = **12,000 fits**.

### Differences from the original script

| | original | here |
|---|---|---|
| Pathfinder start | one run, `effective_N` initialised at the truth | three dispersed starts, no truth |
| ELBO | best-iteration ELBO parsed from Stan's stdout | mean log importance weight over pooled draws |
| Model constants | as written | normalised (adds the `times ≥ 1` truncation; Pathfinder already kept `~` constants — measured, < 0.5 nats) |
| IBD likelihood | normal on IBD fractions with a per-haplotype jackknife SE (Poisson zero term for empty bins) | Poisson on segment counts in every bin (`*_Nfixed_poisson`) |
| Final IBD epoch | `t_end = T_max` | `t_end = t_root + T_max`, as in the Nsmooth models (no-op at these depths) |
| Blocks | 50 × 50 cM, subsamples drawn sequentially | 100 × 25 cM, per-length seed keys |
| IBD sources | two separate scripts and simulations | the same blocks, both sources |

## Outputs (`summary/`)

- `report.md` has, per IBD source, the median ΔELBO and the share of subsamples
  with ΔELBO > 0 for each contrast, plus the share where T_true beats **both**,
  by model and length.
- `contrasts.csv` (per subsample), `fits.csv` (per fit: ELBO, logZ, ESS, k,
  time) and `parameter_estimates.csv` (every mapped parameter of every fit,
  with truth and relative error).
- Figures, one per IBD source: `win_rate_*` (the headline), `delbo_*`
  (distributions), and `parameters_*` (T_true relative error for all 7 times,
  2 fractions and Ne).  `topology_candidates` draws T_true, T_alt1 and T_alt2 with every
  event, fraction and branch labelled by the same names.

## Running

```sh
python execute_full.py --jobs 8                     # compile, simulate, fit, summarise
python execute_full.py --jobs 8 --subsample-limit 1 # pilot: one subsample per length
python -B -m unittest test_design -v
```

Fits run in replicate checkpoints (1, 5, 20, 50, 100 subsamples per length)
and the summary is rebuilt at each one.  Every stage resumes from saved files,
and a manifest guard refuses to mix code or config versions.  Manual stages
are `plan`, `compile`, `simulate [--block-slice k/n]`, `fit [--task-slice k/n]
[--replicate-limit N] [--only-cm ...]`, `summarize` and `calibrate` (NUTS on
T_true at `--calibration-cm`).

## Caveats

Subsamples overlap (up to 40 of the 100 blocks), so shares measure sensitivity
to genome selection, not independent-replicate rates.  The Nfixed models
assume one Ne everywhere, which is true here, so this study tests
identifiability of the two events, not robustness to Ne misspecification (see
`growth_asymmetric_loop`).  Pathfinder is mode-seeking.  Fits that fail the
ESS or Pareto-k screen are kept and counted, never dropped.
