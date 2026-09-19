# Validation (2026-09-10)

- Five scientific-invariant tests passed: 20 genome subsamples with 15 diploid
  individuals per population, correct 2/15/20/40/100/300 event sequence,
  29 orders of 22 graph shapes, separate loop/external fraction indexing, and
  matching semantic parameter/zero-residual calculations.
- Both Poisson models compiled with the explicit two-admixture candidate.
- Small integration test: three 10 cM blocks per scenario, two diploid
  individuals per population, one 20 cM subsample; actual msprime and hap-IBD.
- All 232 candidate fits completed (29 orders × 2 variants × 2 sources ×
  2 scenarios), with compressed posterior draw archives for every fit.
- Parameter/residual tables and 18 PNG/PDF figure pairs were generated.
  Residual figure layout was visually inspected.

Smoke-test results are separate temporary artifacts, not outcomes of the
configured full 20 × 1500 cM study. Full production progress is recorded in
`runs/default/status.json`; no biological conclusion is established by these
execution checks.
