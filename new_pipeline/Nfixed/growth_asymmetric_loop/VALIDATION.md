# Validation (2026-09-19)

Machine: macOS, 12 cores, `genetics_env` (msprime 1.3.4, tskit 0.6.4,
cmdstanpy 1.3.0 + CmdStan 2.37.0, arviz 0.23.4 installed for this study),
Java 21 and hap-ibd 1.0.rev20May22 from the same env, Apple clang 15.

## Design invariants (`test_design.py`, 5 tests, pass)

- **Generating sizes in msprime**, read back from `DemographyDebugger` for both
  scenarios: a = 150,000 at t=0 and 15,000 at t=50⁻; b = 60,000 at t=0 and
  15,000 at t=10⁻; loop1 / loop2 = 15,000 / 5,000 mid-loop; b1 / b2 =
  15,000 / 5,000; root 15,000. Event times 10, 25, 30, 50, 100, 300. The
  `constant` control reads 15,000 everywhere.
- **Composite effective size of the loop**: 1 / Σ p_k²/N_k = 20,000 under
  growth (0.75²/15000 + 0.25²/5000) and 24,000 under constant — the loop
  raises the effective size even with equal branch sizes, because lineages in
  different branches cannot coalesce.
- **Harmonic mean of an exponential branch** matches N₀ ln G / (1 − 1/G)
  (38,376 for a under growth) to 0.1%; the constant control gives 15,000 for
  every reference statistic.
- 29 candidate orders over 22 graphs, 2 explicit-loop orders, 4 correct
  candidates; every node of a correct candidate maps to a generating or
  backbone branch; candidate slices partition the set exactly.
- Residual-tilt slope: 0 on flat residuals, 4/18 per cM on a residual rising
  4 units across the 18 cM bin range, correlation 1.
- `summaries()` carries the per-node references and the log Ne_IBD/Ne_SNP
  contrast (log 2 for a 2:1 input).

## Models and toolchain

- Both normalised Poisson models compile in this env (15 s each); the
  normaliser found exactly 7 and 11 sampling statements as expected.
- msprime growth mechanism verified independently before the study was
  written: leaf `growth_rate` and ancestral `population_parameters_change` at
  the young end reproduce 150,000→15,000 over [0,50] and 30,000→15,000 over
  [50,100] exactly.

## End-to-end smoke test (scratch config)

Three 10 cM blocks, 2 diploids per population, one 20 cM subsample, real
msprime and real hap-IBD; graphs 1, 15, 22 fitted with both models and both
IBD sources for both scenarios (40 fits). All stages completed, tables were
written and all 20 figure families rendered; layouts were inspected visually
(ne_trajectory, spectrum_fit, evidence_benefit, ne_ibd_vs_snp, residual_tilt,
component_likelihood, topology_generating). One bug found and fixed: the
spectrum panel's y-limit collapsed to zero on a pair with no segments.
Smoke-test values are not results.

## Production pilot

Sequential simulation measured at ~92 s per 50 cM block, i.e. ~2.6 h for the
100 blocks before any fit; the simulate stage was parallelised over blocks
(`--block-slice`) and the MRCA scanner precompiled in `compile` so workers
never race on it. Pilot (first subsample, all 232 candidate fits) results are
appended below when complete.

### Pilot result (first subsample, 232 fits, 2026-09-20)

- Simulation: 100 blocks in ~24 min on 8 workers (was ~2.6 h sequential).
  Two parallel-worker races found and fixed on the way: `save_json` used one
  shared temp path (per-PID temp now), and arviz's own import writes a daily
  stamp through a shared temp path (the supervisor imports arviz once before
  launching workers so they never write it).
- Fits: 232/232 completed, 0 failed, median 31 s (shared) / 33 s (separate).
- **0 of 232 pass the importance-sampling screen**: pooled ESS 1–21 of 12,000
  draws, Pareto k 1.5–4 on every fit. With k > 1 the weights have infinite
  mean and the IS logZ is one draw's weight. ELBO (mean log-weight) is
  therefore the primary ranking statistic in `visualize.py`; logZ is kept as
  a secondary column. Same degeneracy Pathfinder showed on the real data.
- One-subsample readout, all conditions: backbone recovered under ELBO;
  ΔELBO(separate − shared) = −45/−24 (growth, true/hap-IBD) and −67/−78
  (constant) on g15; g22 under `constant`+separate ranks first, under
  `growth` second/fourth. log(Ne_IBD/Ne_SNP) on `a`: +0.81 under growth vs
  +0.10 under constant. No conclusion from n = 1; full run launched.

### Full run (10 subsamples, 2,320 fits, 2026-09-20)

- 2,320/2,320 completed, 0 failed, 0 pass the ESS/Pareto-k screen (as in the
  pilot); 165 min wall-clock on 8 workers. All 80 topology comparisons complete.
- Topology (ELBO): shared Ne selects the backbone g15 in 40/40 comparisons and
  never the explicit loop (ELBO(g22) − ELBO(g15) = −12 to −15 nats in every
  subsample). Separate Ne selects the explicit loop g22 — the true graph — in
  70%/50% (growth, true/hap-IBD) and 90%/100% (constant) of subsamples, with
  within-model ELBO(g22) − ELBO(g15) = +9.5/+17.0 and +33.7/+21.1; its backbone
  recovery is 90%/70% (growth) and 90%/100% (constant), i.e. 4 losses, 3 of
  them with hap-IBD under growth.
- Paired on the same graph and data, separate − shared ELBO is negative in
  80/80 comparisons: −57/−69 (growth, g15), −93/−75 (constant, g15), −33 to
  −47 (g22).
- Branch a under growth: shared 39,809 vs harmonic mean 38,376 (arithmetic
  58,630); separate IBD 40,121, separate SNP 17,192. log(Ne_IBD/Ne_SNP) on a:
  +0.83 [+0.76, +0.89] growth vs −0.09 [−0.14, −0.03] constant (true IBD);
  +0.86 vs −0.05 with hap-IBD. On collapsed b: ≈ 0 in both scenarios.
- Small branches: shared recovers b.2 = 7,559 (true 5,000) and right = 9,691
  (true 15,000); separate IBD gives 25,750 / 24,675 for the same branches.
- Event times on g15 (median, growth/true): shared 35 / 56 / 131 / 260 for
  b_split / left / right / root (truth 30 / 50 / 100 / 300); separate 28 / 66 /
  64 / 339 — `right` collapses onto `left`.
- Residual tilt is weak evidence: a–a −0.018 (growth) vs 0.000 (constant), but
  c–c, which has no growth, differs by a similar amount between the two pools.
