# b-loop then admixture, with growth and asymmetric branch sizes

Third study in this series, after `hidden_loop_robustness` (hidden loop in c)
and `b_loop_then_admixture` (recent loop in b, constant Ne everywhere).  Same
question — does a **separate Ne for the SNP and IBD components** buy anything
when the fitted model is misspecified? — but now the misspecification is in
**Ne itself**, not only in the omitted loop.  Every candidate has a constant Ne
per branch; the truth has exponential growth and 3:1 size asymmetries.

## Generating history

Populations `a, b, c`, 15 diploid (30 haploid) samples each.  Backwards in time:

| Generations | Event |
|---:|---|
| 10 | b splits into loop1 (0.75) and loop2 (0.25) |
| 25 | loop1 and loop2 merge into anc_b |
| 30 | anc_b splits into b.1 (0.75, toward a) and b.2 (0.25, toward c) |
| 50 | a and b.1 merge |
| 100 | b.2 and c merge |
| 300 | root |

The loop is at 10–25 rather than the earlier study's 2–15 so that b's leaf
branch is long enough for growth to mean anything.

Haploid Ne, scenario **`growth`** (the experiment):

| branch | Ne (young end → old end) |
|---|---|
| a | 150,000 → 15,000, exponential across [0, 50] (10×) |
| b (leaf) | 60,000 → 15,000, exponential across [0, 10] (4×) |
| loop1 / loop2 | 15,000 / 5,000 (3 : 1) |
| anc_b | 15,000 |
| b.1 / b.2 | 15,000 / 5,000 (3 : 1) |
| c, left, right, root | 15,000 |

Scenario **`constant`** (the control): identical topology and times, every
branch 15,000.  Both scenarios contain the loop; the control isolates what the
Ne misspecification alone does.  All values live in `config.json`
(`ne.growth_fold_a`, `ne.growth_fold_b`, `ne.ratio_loop`, `ne.ratio_split`).

Growth is implemented in a local msprime builder (`study.msprime_demography`)
so no shared code changes and the other studies' manifests stay valid.  For an
ancestral branch on [t0, t1] a `population_parameters_change` at t0 anchors the
young-end size and rate; growth is switched off at t1.  Sizes at every branch
end are checked against msprime's `DemographyDebugger` in `test_design.py`.

## What is fitted

Identical search to `b_loop_then_admixture`: the 21 loop-omitting three-leaf
graphs (3 trees + 18 one-admixture shapes, all 27 valid event orders) plus the
explicit generating loop graph 22 (2 orders) — **29 ordered candidates**, each
fitted with `mixed_model_Nsmooth_poisson` (shared Ne) and
`mixed_model_Nsmooth_poisson_separate_ne`, to true IBD and to hap-IBD, for both
scenarios.  Graph 15 `((a,b.1),(b.2,c))` is the correct collapsed backbone.
Three Pathfinder starts per fit (4 paths × 1,000 draws).  Stan copies in the
run's `models/` restore every dropped constant so evidence is comparable across
graph sizes and across shared/separate.

**Graphs are ranked by ELBO** (mean log importance weight, averaged over event
orders), not by the importance-sampling logZ the earlier studies used.  The
pilot showed why: pooled importance ESS of 1–21 out of 12,000 draws and Pareto
k between 1.5 and 4 on every one of 232 fits.  With k > 1 the weights have
infinite mean, so the IS logZ is effectively the single largest weight.  The
ELBO is a low-variance lower bound and is what the study compares; logZ, ESS
and k are still recorded per fit and per start, and `reliable` (ESS ≥ 100,
k ≤ 0.7) is reported honestly — expect it to be false almost everywhere.

## Statistics beyond topology recovery

Under `growth`, times and sizes are **not** recovery targets — the truth is
outside every candidate's model class.  The analysis asks instead how each
model absorbs the misspecification and whether the second Ne trajectory helps:

- **`evidence_benefit`** — paired logZ and ELBO (separate − shared) on the
  same subsample and the same graph, for g15 and g22.  Constants are
  normalised, so a positive value means the second trajectory is supported
  after paying for its parameters.  `paired_ne_comparison.csv` also records
  rescue / loss of the correct backbone.
- **`ne_trajectory`** — the true effective Ne(t) of each backbone branch with
  the fitted constants over it.  The collapsed b branch is a composite: leaf
  growth, then the loop's `1 / Σ p_k² / N_k` (two lineages coalesce only if
  they fall in the same branch; the same quantity governs drift), then anc_b.
  The harmonic mean is what accumulated drift implies (duration / harmonic),
  so it is the natural SNP-side reference; the arithmetic mean is also drawn.
- **`ne_ibd_vs_snp`** — `log(Ne_IBD / Ne_SNP)` per branch from the separate
  model.  Growth should push this above zero on a and b: IBD is dominated by
  the recent, large end, SNP drift by the harmonic mean, which the small old
  end controls.  The constant control should sit at zero.
- **`residual_tilt`** — least-squares slope of the Pearson residual against
  segment length, per pair.  A tilt is the fingerprint a constant Ne leaves
  when the truth changed size inside the branch.
- **`spectrum_fit`** — observed counts (one subsample, exact Poisson 68%
  bars, zero-count upper limits) against each model's fitted expected counts,
  with median ± IQR Pearson residuals across all subsamples below.
- **`component_likelihood`** — lp and chi²/n by data type.
- **`parameters_backbone`** — event times and fraction on g15, descriptive.

`summary/report.md` tabulates recovery, the paired benefit, fitted-vs-true
branch sizes (young / old / harmonic / arithmetic vs shared / IBD / SNP) and
the tilt.  `ne_estimates.csv` carries every fitted size with its reference.

## Running

Requires `genetics_env` (msprime, tskit, cmdstanpy, arviz, scipy, matplotlib),
CmdStan, a C++ compiler, and hap-IBD with Java; the last two ship inside the
conda env and `execute_full.py` defaults to them (`JAVA`, `HAPIBD_JAR` override).

```sh
python execute_full.py --jobs 8                 # compile, simulate, fit, summarise
python execute_full.py --jobs 8 --subsample-limit 1   # pilot
```

Fits run as parallel workers over disjoint candidate slices; every stage
resumes from saved files.  Manual stages: `plan`, `compile`, `simulate`,
`fit [--candidate-slice k/n]`, `summarize`, `calibrate`.  Run folders are
ignored by Git; a manifest guard refuses to mix code/config versions.

```sh
python -B -m unittest discover -s . -p 'test_design.py' -v
```

## Caveats that travel with the results

Subsamples overlap (30 of 50 blocks), so rates are sensitivities to genome
selection, not independent-replicate frequencies.  Rankings are composite-
likelihood importance estimates, not calibrated probabilities.  Pathfinder
is mode-seeking; a completed fit that fails the ESS / Pareto-k checks is kept,
labelled unreliable, and never silently dropped.  Sizes on collapsed branches
are interpretive: the fitted constant is being asked to stand in for a curve.
