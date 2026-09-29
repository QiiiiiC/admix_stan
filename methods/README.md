# methods

Shared, plotting-free machinery for admixture-graph simulation studies that fit
IBD segments and SNP allele-frequency covariance with Stan.  Extracted from
`new_pipeline/` and the `growth_asymmetric_loop` study; every study under
`simulation_study/` imports it.

```python
import sys; sys.path.insert(0, "<repo root>")
from methods import ibd, snp, blocks, models, fitting, stan_data, simulate
```

Requires `genetics_env` (msprime, tskit, cmdstanpy + CmdStan, arviz, scipy,
networkx for `graph_search`), a C++ compiler, and Java + hap-IBD for called IBD.

## Pipeline

| Step | Module | Main entry points |
|---|---|---|
| Describe a history | `demography` | `DemographicTopology` (`add_merge_event`, `add_admixture_event`, `set_*`), `event_times`, `admixture_fractions`, `replay`, `valid_event_orders`, `clades` |
| Candidate graphs | `topologies` | `enumerate_all` (21 three-leaf graphs), `ordered_candidates` (every valid event order) |
| | `graph_search` | `enumerate_topologies(leaves, max_admix)` for n leaves |
| Simulate | `simulate` | `to_msprime(dem, sizes={node: (young, old)})` with exponential growth, `add_branch`, `simulate_block`, `block_seed` |
| IBD | `ibd` | `true_ibd` (strict MRCA spans, C++ `mrca_scan.cpp`), `hapibd`, `bin_edges`, `pair_counts` |
| SNP | `snp` | `snp_summaries` (per sub-block sums), `snp_covariance` (ratio-of-sums W, jackknife SE) |
| Observed data | `blocks` | `selections` (blocks without replacement; extra key e.g. genome length), `load_blocks`, `aggregate` |
| Stan data | `stan_data` | `stan_data(dem, obs, spec, edges, n_haploid, T_max)`; `topology_data`, `smooth_data` |
| Models | `models` | `MODELS` registry, `normalized_model(spec, folder)` |
| Fit | `fitting` | `pathfinder` (multi-start, pooled importance weights), `nuts`, `initial`, `parameter_summary` |
| Check fit | `diagnostics` | `fit_summary`, `residual_tilt`, `trajectory_reference`, `true_parameters` |
| Run at scale | `parallel` | `Supervisor.stage(...)` over `--block-slice`/`--candidate-slice` workers, `in_slice` |
| Housekeeping | `utils` | `save_json` (atomic, per-PID), `open_run` (manifest guard), `file_hashes`, `quiet` |

## Models (`stan/`)

| Registry name | Data | Ne | IBD likelihood |
|---|---|---|---|
| `ibd_Nfixed`, `snp_Nfixed`, `mixed_Nfixed` | IBD / SNP / both | one `effective_N`, gamma(4, 4/15000) | normal on `ibd_hat`/`ibd_se`, Poisson zero term for empty bins |
| `ibd_Nsmooth`, `snp_Nsmooth` | IBD / SNP | per-branch log-normal random walk | normal (as above) |
| `mixed_Nsmooth_poisson` | both | random walk | Poisson on raw counts, every bin |
| `mixed_Nsmooth_poisson_separate` | both | independent IBD and SNP random walks | Poisson |

`normalized_model` compiles a copy whose `lp__` keeps every constant: `~`
statements become explicit `*_lpdf`, and truncation constants (`times >= 1`
under an exponential prior, half-Normal scales) are derived from the source and
added.  Posteriors are unchanged.  This is what makes ELBO/logZ comparable
between graphs with different numbers of events.

## Conventions

- Sizes: `DemographicTopology` stores **diploid** `node.ne`; Stan and
  `simulate` use **haploid** sizes (true haploid = `2 * node.ne`).
- Stan node index = insertion order of `dem.nodes`; event k = `cumulative_times[k]`;
  `admixture_fractions[i]` = fraction from the FIRST parent of the i-th admixture.
- IBD arrays are `[bin, i, j]`, symmetric, populations in the order passed;
  exposure is haplotype pairs (`n(n-1)/2` within, `n_i n_j` between).
- `blocks.aggregate` output (`ibd_count`, `ibd_hat`, `ibd_se`, `w_hat`,
  `w_se`, `cm`) serves every model; `ibd_se` is the SE across the chosen
  blocks, so a subsample needs at least two blocks.
- Pathfinder evidence: rank by ELBO; logZ, ESS and Pareto k are recorded but
  k > 1 is common, in which case logZ is one draw's weight.

## Tests

```sh
python -B -m unittest methods.test_methods -v     # from the repository root
```

Covers the DSL truth accessors, event orders, the three-leaf search, msprime
growth against `DemographyDebugger`, SNP covariance, block aggregation, the
normaliser on every model, Stan data keys against every model's `data` block,
and the harmonic-mean reference.
