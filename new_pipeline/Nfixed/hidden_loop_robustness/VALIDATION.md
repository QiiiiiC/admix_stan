# Implementation validation

Validated locally on 2026-09-08 with the existing `stan_env` environment:
NumPy 2.3.3, msprime 1.3.4, tskit 0.6.4, CmdStanPy 1.3.0 and CmdStan 2.37.0.
The real hap-IBD JAR in `/Users/qi/software/hap-ibd/` ran using the locally
cached OpenJDK 25.0.2 executable.

- Eight unit tests passed, including strict MRCA-span equivalence, binning,
  raw counts and exposure, both event orders, SNP jackknife, and HBD inclusion.
- All four normalized non-grid Stan model copies compiled.
- A deliberately small integration experiment generated three 10 cM blocks
  per scenario with two diploid samples per population, then selected two
  blocks (20 cM). It used one Pathfinder start/path and 100 requested draws.
- Both true IBD and actual hap-IBD were processed for both scenarios.
- All **432 candidate fits** completed: 27 event orders × four models × two
  IBD sources × two scenarios. Reports generated all **16 complete rankings**.
- The report also correctly marked an earlier partial-graph run incomplete.
- Compressed importance weights remained available after successful-fit
  Pathfinder CSV cleanup; zero such CSVs remained in the completed smoke run.
- The optional four-chain MCMC stage ran with 30 warmup and 30 sampling
  iterations. Its maximum R-hat was 1.08811 and it was correctly flagged
  unreliable. This checks execution and diagnostics, not interval calibration.
- At identical parameter values, the original and normalized shared-Poisson
  model log densities differed by the analytically expected constant
  (-2.18700959843752 versus -2.1870095984375117). Their gradients agreed.

These checks establish that the implementation executes; they do **not**
identify a winning model or establish statistical calibration. The configured
full 1500 cM, 100-subsample study has not been run. Its settings remain in
`config.json`; smoke settings and results were kept separately in temporary
directories and did not alter that configuration.
