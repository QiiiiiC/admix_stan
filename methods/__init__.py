"""Shared methods for admixture-graph simulation studies (IBD + SNP, Stan).

Modules, in pipeline order:

  demography    DemographicTopology DSL; replay / valid event orders / clades
  topologies    exhaustive three-leaf search (21 graphs, every event order)
  graph_search  general n-leaf, k-admixture enumeration (networkx)
  simulate      DSL -> msprime (with per-branch exponential growth); one block
  ibd           true IBD (MRCA spans, C++) and hap-IBD calls, binned by length
  snp           per-block SNP covariance sums; ratio-of-sums W with jackknife SE
  blocks        block subsampling and aggregation into observed data
  stan_data     Stan data for IBD-only / SNP-only / mixed, Nfixed / Nsmooth
  models        Stan model registry; normalised copies for evidence comparison
  fitting       multi-start Pathfinder + importance sampling; NUTS
  diagnostics   fit summaries, residuals, residual tilt, true-size references
  parallel      worker supervisor for sliced, resumable stages
  utils         atomic JSON, hashing, run-manifest guard

No plotting lives here.  Import submodules directly, e.g.
`from methods.ibd import true_ibd`.
"""
