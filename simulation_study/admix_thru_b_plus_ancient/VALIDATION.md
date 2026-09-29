# Validation (2026-09-28)

Machine: macOS, 12 cores, `genetics_env` (msprime, tskit, cmdstanpy + CmdStan,
arviz), Java + hap-ibd 1.0.rev20May22 from the same env.

## Design invariants (`test_design.py`, 5 tests, pass)

- msprime history: event times 20/60/100/500/700/900/1000, every population at
  haploid 10,000 with no growth, admixture proportions 0.5/0.5 on a and 0.7/0.3
  on c with cP1 first.
- Candidates: T_true = the generating graph's events; 7/5/5 events and 2/1/1
  admixtures; all valid.
- Semantic mapping: T_true's 7 times and 2 fractions equal
  `true_topology().event_times()` / `admixture_fractions()` with correct
  indices; T_alt1 maps only loop_close, the ancient event and the merges above
  it; T_alt2 only the loop.
- Task list: 12,000 unique fits; SNP-only only under source `none`; 8 slices
  partition it; a 1-subsample pilot is 120 fits.  Block selections have the
  right size, are deterministic and differ between lengths.
- Stan data for every candidate × model; summaries recover zero error on
  truth-valued draws.

`methods` checks (`python -B -m unittest methods.test_methods`, 9 tests,
pass) include Stan data keys equal to each Nfixed model's `data` block and the
normaliser finding 3 / 4 / 4 sampling statements.

## Engine constant check

The original study ranked by the ELBO Pathfinder prints.  Since T_true has 7
events and the alternatives 5, dropped `log(0.01)` constants of the
exponential time priors would have biased ΔELBO by ~9.2 nats.  Measured on
snp_Nfixed: raw vs normalised model differ by 0.41 nats (stdout ELBO) and 0.32
(lp__), against 76 expected if `~` constants were dropped — Pathfinder keeps
them.  The original numbers are therefore not biased by constants; only the
`times ≥ 1` truncation (0.01 per event) is missing there.

## End-to-end smoke test (scratch config)

4 × 10 cM blocks, 4 diploids per population, lengths 20 and 30 cM, 2 subsamples,
2 starts × 200 draws, real msprime and real hap-IBD: 60/60 fits ok; summary
tables, report and all 6 figures written and inspected (end-label collisions
and share labels on whiskers fixed).  Already at 20–30 cM the division of labour
appears: IBD-only wins the recent contrast, SNP-only and Mixed the ancient one.
Smoke-test values are not results.

## Pilot 1 (normal IBD likelihood) — stopped, superseded

The first production launch used `ibd_Nfixed` / `mixed_Nfixed`.  Their normal
likelihood needs an IBD SE, which was taken as the SE across the chosen blocks.
The first checkpoint (120 fits, 0 failed, median 10 / 5 / 15 s for IBD / SNP /
mixed) showed IBD-only ΔELBO(T_true − T_alt1) = +1011 at 50 cM (true IBD) and
−91 (hap-IBD), against +11 at 100 cM.  Cause: with 2 blocks the SE has one
degree of freedom; 7 cells at 50 cM had SE/mean between 1.5% and 10% of the
Poisson counting error 1/sqrt(k) (worst at 100 cM: 0.23×; at 1000 cM: 0.80–
1.49×).  The run was stopped and kept as `runs/pilot_normal_likelihood`.

Fix: `ibd_model_Nfixed_poisson.stan` and `mixed_model_Nfixed_poisson.stan`
(in `methods/stan/`).  These are the Nfixed models with the Poisson count
likelihood of the Nsmooth Poisson models, `ibd_count` in place of
`ibd_hat`/`ibd_se`, component `lp_*`/`chi2_*` generated quantities, and the
root-anchored final epoch.  `methods.stan_data` now passes exactly the
variables a model's `data {}` block declares.  The rerun smoke test (60/60 ok)
gives IBD-only ΔELBO ≈ 0 at 20–30 cM with 4 diploids, as the data warrant.

## Full run (Poisson IBD; 100 subsamples × 8 lengths, 12,000 fits, finished 2026-09-28)

12,000/12,000 completed, 0 failed; 748 pass the ESS/Pareto-k screen (ELBO is
the ranking statistic, as designed).  `summary/` is copied out of the
git-ignored `runs/default/summary`.

- Recent contrast (T_true − T_alt1), true IBD: IBD-only median ΔELBO +2.8 →
  +134 over 50 → 1000 cM (wins 84% → 100%); Mixed +4.2 → +149 (83% → 100%);
  SNP-only negative throughout (median −4.7 → −1.4; wins 1–34%).
- Ancient contrast (T_true − T_alt2): SNP-only and Mixed +75 → +1370, 100% at
  every length.  IBD-only sits at a flat +0.6 to +0.9 nats that does NOT grow
  with genome length (wins 72–89%): a fixed offset between the two graphs, not
  detection — real signal scales with length, as every other panel does.
- T_true beats both alternatives: Mixed 83% (50 cM), 99% (100 cM), 100% from
  150 cM; hap-IBD 89% / 99% / 100%.  IBD-only's "both" (68–88%) inherits the flat
  offset above; SNP-only 1–34%.
- True IBD and hap-IBD agree to within a few nats on every contrast.
- Recovery (Mixed, T_true): Ne within ~2%, root/right/left within ~5% at
  ≥ 500 cM; the ancient fraction is +8% biased under both SNP-only and Mixed,
  and IBD-only puts it at −38% (IBD cannot see it).  SNP-only mis-sizes the loop
  (loop fraction −55%, donor merge −28%), which Mixed corrects to IBD's values.

## Why SNP-only's Ne is off (checked 2026-09-28)

All three models are Nfixed: one free `effective_N` for every branch, never
given the truth.  SNP drift on a branch is duration / Ne, so SNP alone identifies
only the ratio.  Medians over the 100 true-graph fits:

| model | cM | Ne / truth | t_root / truth | (t_root/Ne) / truth | (t_left/Ne) / truth |
|---|---:|---:|---:|---:|---:|
| SNP-only | 50 | 0.870 | 0.866 | 1.007 | 1.007 |
| SNP-only | 200 | 0.869 | 0.865 | 1.000 | 1.009 |
| SNP-only | 1000 | 0.913 | 0.920 | 0.998 | 1.002 |
| IBD-only | 1000 | 0.989 | 0.928 | 0.936 | 0.994 |
| Mixed | 1000 | 0.996 | 1.000 | 1.003 | 1.037 |

SNP-only gets the drift exactly and slides along the (t, Ne) ridge, so Ne and
the times are low together.  IBD supplies the scale; Mixed recovers both.
