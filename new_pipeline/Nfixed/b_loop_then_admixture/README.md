# Recent b loop followed by admixture

New experiment, independent of `hidden_loop_robustness` (the earlier c-loop study).
The user's **20 samples means 20 simulation subsamples**, not 20 individuals.
Each population has **15 diploid individuals (30 haplotypes)**. All generating branches have haploid
Ne 15,000, passed to msprime as diploid size 7,500.

Backward in time:

| Time (generations) | Event |
|---:|---|
| 2 | b splits into loop1 and loop2, fraction 0.75 toward loop1 |
| 15 | loop1 and loop2 merge into anc_b |
| 20 | anc_b splits into b.1 and b.2, fraction 0.75 toward b.1 |
| 40 | a and b.1 merge |
| 100 | c and b.2 merge |
| 300 | Both ancestors merge into root |

The no-loop control deletes events 2 and 15: b itself admixes at 20.
There is no loop in c. Fraction 0.75 for the loop and the later merges at
40/100/300 are retained defaults; both fractions are inferred when fitting
the explicit-loop graph, not fixed at their generating values.

## Data and search

For each scenario, simulate **50 independent 50 cM blocks** once. Each of
20 subsamples selects **30 distinct blocks = 1500 cM**, without replacement.
The same block selection supplies SNPs and both IBD sources to all fits.
Selections overlap across subsamples; no independent-replicate confidence
claims are made. No held-out blocks. Both scenarios use their own simulated pool.

Uniform 0.5 cM bins span `(2,2.5], ..., (20,20.5]`. True IBD uses strict
continuous MRCA-node spans. hap-IBD uses simulated phased genotypes and reads
both IBD and HBD output to match all 435 within-population and 900 between-
population haploid pairs. The SNP ascertainment, drift normalization and
physical-block jackknife match the earlier study: ancestral SNPs (mutation
node at least as old as the root), global MAF > 0.05, 10 cM SE sub-blocks.
This ancestral ascertainment uses simulation knowledge; it is not a real-data
filter. Caller input is unfiltered phased SNP data.

Fit **Poisson shared Ne and Poisson separate Ne**, both non-grid Nsmooth
models from `real_data`, to true IBD / hap-IBD and b-loop / no-loop scenarios.

- Graphs 1–21: the entire existing loop-omitting search, including three trees
  and 18 one-admixture shapes, with all 27 valid chronological orders.
- Graph 15 is `((a,b.1),(b.2,c))`, the correct collapsed backbone.
- Graph 22 adds the explicit generating b loop. It has two orders of the
  independent a/b.1 and c/b.2 merges, giving **29 ordered candidates total**.
- Graph 22 is fitted even to the no-loop control to assess spurious loop selection.
- Total: **29 × 2 models × 2 IBD sources × 2 scenarios × 20 subsamples = 4,640 fits**.

Each fit uses three Pathfinder starts, each with four paths and 1,000 requested
draws per path. Graph scores average evidence across orders with equal
conditional weights, using the same fully normalized Stan copies as the earlier
study. The original Stan sources remain unchanged. Topology weights are based
on the existing composite likelihood, not calibrated probabilities of biological
truth. Evidence and parameter approximations retain ESS/Pareto-k diagnostics.
Adding separate Ne or an unobserved loop can increase nonidentifiability;
successfully completing a fit is not equivalent to passing these diagnostics.

## Running and progress

The active production run was capped at **10 subsamples (2,320 fits)** on
2026-09-10 at the user's request. The immutable manifest retains the original
20 selections; only the first 10 are fitted. Resume this capped run with:

```sh
python execute_full.py --subsample-limit 10
```

From this folder, with the existing `stan_env` Python:

```sh
python run_study.py plan
python run_study.py simulate --hapibd-command 'java -jar /path/to/hap-ibd.jar'
python run_study.py fit
python run_study.py summarize
```

`execute_full.py` uses the locally available Java/JAR paths, runs the complete
workflow, and writes `runs/default/status.json`, `simulate.log`, `fit.log`, and
`summarize.log`. It fits the first subsample across all conditions, renders
figures, then extends to the second and so on through 20. Existing fits and
simulation blocks are reused. Figures therefore first appear **after simulation
and the first subsample's fits**, not immediately upon launch.

Every invocation validates configuration and source hashes against the run
manifest. Each new experiment uses a new `--out` directory. Tests:

```sh
python -B -m unittest discover -s . -p 'test_design.py' -v
```

## Saved information and visual output

Every candidate's `result.json` includes topology and event identities, per-start
evidence/ELBO/ESS/Pareto-k, runtime/failures, all event-time/fraction/Ne intervals,
semantic backbone and loop parameters when comparable, observed summaries,
prediction means/intervals, component likelihood summaries, and residual arrays.
`posterior_draws.npz` retains all sampled parameters, transformed sizes/times,
component scores and raw importance weights for later analyses. Full prediction
draws can be regenerated from these parameters; their intervals are already
saved. Large Stan CSVs are removed after saved summaries/draw archives unless
`keep_stan_csv` is enabled. Tree sequences, caller inputs/logs, block sufficient
statistics, random seeds and selected blocks are retained.

The two admixture fractions are distinguished explicitly: in graph 22,
`admixture_fractions[0]` is the loop and `[1]` is the external b admixture.
The loop fraction is plotted as `max(f,1-f)` because its source labels exchange.
The later fraction is always oriented toward a. Loop times/fractions have
**no true value** in the no-loop control; their truth fields are null there.

`summary/` saves CSV/JSON tables for later use:

- `search_rankings.csv` and `search_selection.csv`: restricted loop-omitted
  search (21 graphs) AND augmented search (22), exact graph vs backbone recovery,
  evidence gaps and ELBO sensitivity. Incomplete comparisons are marked.
- `parameter_estimates.csv`, `parameter_metrics.csv`: raw and reliable-only
  estimates/errors/coverage, with explicit candidate/event/component identities.
- `residuals.csv`: raw fraction/count, Pearson/deviance IBD, raw/standardized SNP
  residuals for every fitted candidate. Zero count bins are retained.
- `fit_diagnostics.csv`, `fit_index.json`, `summary_status.json`: execution and
  numerical diagnostics, plus paths to each fit and draw archive.

Only **parameter and residual figures** are rendered, in `summary/figures/`,
as PNG and PDF. No topology ranking plots are produced.

- `parameters_matched_backbone_*`: compare graph 15 with graph 22 at the
  best-evidence event order within each graph. These are conditional graph comparisons.
- `parameters_selected_topology_*`: best graph from each search. If a winning
  graph lacks a directly matching backbone parameter, it is omitted and its
  smaller sample count is displayed rather than mislabeling an event.
- `loop_parameters_*`: explicit-loop opening, closure and folded loop fraction.
- `Ne_*`: shared and separate IBD/SNP branch-size estimates for each graph.
- `residuals_*`: six IBD count-residual spectra plus SNP residuals, contrasting
  omitted/explicit loop and shared/separate Ne. Curves are subsample medians;
  bands/bars are IQR, not posterior or independent-replicate confidence intervals.

All parameter boxes show variation of weighted posterior means across genome
subsamples, with point counts and generating truth lines where defined.
Unreliable fits remain visible and labeled; plots are descriptive.
