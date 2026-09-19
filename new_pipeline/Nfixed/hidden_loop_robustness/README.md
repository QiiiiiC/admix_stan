# Hidden c-loop: robustness of inferred backbone topology

This study compares the four **non-grid Nsmooth mixed models** in
`new_pipeline/real_data`: Poisson/Normal × shared/separate SNP–IBD Ne.
`Nfixed` here describes the simulated demography: all generating branches
have haploid Ne = 15,000 (msprime diploid size = 7,500). The fitted models
still infer their branch sizes. No existing model or simulation script is modified.

## Design

Present-day populations are `a, b, c`, with 15 diploid / 30 haploid samples
each. Events run backwards from the present:

| Generations ago | Hidden-loop generating event |
|---:|---|
| 12 | c splits into c1 and c2, fraction 0.75 toward c1 |
| 16 | b splits into b1 and b2, fraction 0.75 toward b1 (the a side) |
| 40 | a and b1 merge into left |
| 70 | c1 and c2 merge into anc_c |
| 100 | b2 and anc_c merge into right |
| 300 | left and right merge into root |

The **no-loop control** removes the events at 12 and 70; c joins b2 at 100.
Both generating scenarios therefore have the same target backbone,
`((a,b.1),(b.2,c))`. The hidden loop is excluded from every candidate fit.

- Generate one independent pool of **50 × 50 cM chromosomes per scenario**.
- Draw **30 distinct chromosomes without replacement** for each 1500 cM
  dataset. Repeat 10 times by default; chromosomes can recur across datasets.
- Identical selected chromosome IDs supply SNPs, true IBD, and hap-IBD to every
  model and candidate topology. IDs are also matched across the two scenarios,
  whose underlying ancestry simulations use different random seeds.
- **No held-out, validation, or test blocks.** All observations from the selected
  30 chromosomes are available for fitting.
- Uniform bins `(2,2.5], (2.5,3], …, (20,20.5]` cM, matching the repository's
  lower-open/upper-closed convention. Segments outside these bins are omitted.
- Recombination rate `1e-8/bp/generation`, mutation rate `1.25e-8`;
  50 cM corresponds to 50 Mb. The default `smc` ancestry model follows the
  existing simulation workflow; it is configurable.

**Subsamples are dependent.** Recovery rates describe sensitivity to genome
selection within the simulated pool. They do not estimate independent-genome
coverage with a binomial standard error. No such error bars or tests are used.
`pools > 1` optionally regenerates independent pools; the default remains one
per scenario to save simulation cost. Pool IDs are retained in every result.

## What is compared

The primary endpoint is **exact backbone graph recovery**, ignoring the order
of independent merge events. Also saved: true graph rank, true-versus-best-wrong
score gap, recovery of b as the admixed leaf, and ELBO-ranking sensitivity.

All 21 existing three-population graph shapes are included: three trees and
18 one-admixture graphs. Every valid temporal ordering is represented, giving
**27 ordered candidate models**. This includes both orders of independent
merges and cases where the two non-admixed populations merge before admixture.
The target is graph 15, with ordered models `g15_o1` and `g15_o2`; only `o1`
has the generating event order. These are graph candidates, not a claim that
the true graph is identifiable from three populations in every dataset.

Within **each model variant and IBD source**, rank graph shapes using the
importance-sampling estimate of marginal likelihood from Pathfinder:

```text
log_weight = lp__ - lp_approx__
logZ(order) = logmeanexp(log_weight across all starts)
logZ(graph) = logmeanexp(logZ(order) across valid orders)
```

Graph shapes receive equal prior weight; orders receive equal conditional
weight within their graph. Averaging instead of maximizing avoids an automatic
advantage for shapes with more permissible orders. ELBOs also combine across
orders with logmeanexp, providing an alternative lower-bound ranking.
All starts are retained rather than selecting the most favorable start.

The generated Stan copies in the run's `models/` directory replace sampling
statements with full `*_lpdf` calls. They also include the normalizers for
`times >= 1` (truncated exponential) and positive sigma/tau (half-Normals).
These changes are **additive constants only**, preserving the original
posterior distributions while making evidence comparable across graph sizes.
Original files in `real_data` remain unchanged.

Pathfinder importance estimates can be unstable. The default checks require
all starts to succeed, no discarded nonfinite weights, pooled importance ESS
at least 100, and Pareto k at most 0.7. Per-start scores and diagnostics are
retained. Recovery is reported both for all complete comparisons and for the
subset where **all 27 candidate fits pass**. Missing fits do not become wrong
topology calls. A zero reliable denominator is shown as unavailable.
These diagnostics do not prove that all posterior modes have been found.
The scores use the existing SNP + IBD **composite likelihood**; they are not
calibrated probabilities that a demographic graph is true.

**Never compare Poisson and Normal logZ/ELBO levels directly:** their IBD
observations are different (counts versus length fractions). Compare their
topology recovery frequencies, ranks, parameter errors and residuals instead.
The `paired_ne_comparison.csv` directly records when separate Ne rescues or
loses a topology call relative to shared Ne on the identical dataset.

Secondary statistics:

- Relative time bias/RMSE for b admixture, left merge, right merge and root;
  absolute bias/RMSE for b's fraction.
- 95% interval coverage and widths, **conditional on fitting the true backbone
  and true order**, not after selecting a topology. Pathfinder intervals are
  approximate. The optional `calibrate` stage uses four-chain MCMC for these
  intervals, with R-hat, ESS, divergence and tree-depth diagnostics.
- IBD fraction residuals by all six population pairs and 37 bins; per-pair
  RMSE; SNP covariance RMSE and standardized residuals. These describe fitted
  data, not independent predictive performance.
- Branchwise log(estimated Ne / 15,000) and posterior log(Ne_IBD / Ne_SNP)
  contrasts. Sizes on collapsed-loop branches are interpretive diagnostics,
  not necessarily directly equivalent demographic parameters.
- Runtime and failure/importance-weight diagnostics for all attempted fits.

## IBD and SNP handling

**True IBD** uses continuous spans with the same MRCA **node**, matching
`ibd_jackknife.calculate_ibd_blocks_mrca`, not ancestry-path IBD from
`ts.ibd_segments`. The small C++ helper scans the same MRCA spans more quickly.
No minimum-span filter is applied before adjacent spans are merged. Tests
compare it against a direct tskit MRCA scan and a path-changing example.
Counts are actual segment counts; they are never reconstructed by dividing
total length by bin midpoint. Entire chromosomes are sampled, so subsampling
does not cut segments or join unrelated chromosomes.

**hap-IBD** receives the simulated phased VCF, including all variants, and a
uniform map. There is no additional phasing/genotyping error layer. Both
`.ibd.gz` **and `.hbd.gz`** are read, because the Stan diagonal exposure
`30*29/2 = 435` includes the two homologues within each diploid individual.
Between-population exposure is 900. Caller settings, logs and JAR checksum
(when using a JAR) are retained. Defaults: min-seed 1 cM, min-output 2 cM,
one thread; other caller parameters retain the installed version's defaults.
See the [hap-IBD output specification](https://github.com/browning-lab/hap-ibd#output-files).

SNP summaries follow the existing simulation's **ancestral-SNP ascertainment**:
retain biallelic, single-mutation sites whose mutation-bearing node is at least
as old as the generating root, with global MAF > 0.05. This uses known
simulation ancestry and is not an implementable real-data ascertainment rule.
Set `snp_ancestral_only: false` to retain qualifying variants of all ages;
that is a separate experimental setting. hap-IBD always receives the unfiltered
variant data.

SNP covariance uses the repository's TreeMix ratio-of-sums and centered
`1/n_haploid` finite-sample correction. There is no LD-pruning step in the
reused block-summary workflow. SEs use a delete-one-**10 cM physical block**
jackknife over the selected chromosomes (150 sub-blocks at 1500 cM), rather
than windows of a fixed SNP count that change physical size with SNP density.
No physical SE block bridges chromosomes. IBD block SEs are descriptive;
the four Stan models use their existing count or theoretical-Normal likelihood.

## Running

Requires Python with `numpy`, `scipy`, `msprime`, `tskit`, `cmdstanpy`, `arviz`,
and `matplotlib`; a working CmdStan installation; a C++ compiler; and hap-IBD
with Java. The existing local `/opt/anaconda3/envs/stan_env/bin/python` has the
Python dependencies. Use your normal Java installation or supply its full path.

Run from this directory:

```sh
python run_study.py plan
python run_study.py simulate --hapibd-command 'java -jar /path/to/hap-ibd.jar'
python run_study.py fit
python run_study.py summarize
```

The complete default design is **4,320 ordered candidate fits**, each with
three multi-path Pathfinder starts. For a first pilot, use the same full pools
and fit only the first subsample; extend later without repeating completed work:

```sh
python run_study.py fit --replicate-limit 1
python run_study.py summarize
python run_study.py fit
```

The pilot report explicitly marks the other 9 subsamples as incomplete.
For MCMC-based parameter intervals on the correct backbone/order:

```sh
python run_study.py calibrate --replicate-limit 1
python run_study.py summarize
```

`calibrate` does not produce evidence estimates or replace topology rankings.
`--warmup`, `--samples`, and `--parallel-chains` control MCMC effort.

All stages resume from saved files. `--retry-failed` retries unsuccessful or
unreliable fits. Tree sequences, compressed per-block sufficient statistics,
hap-IBD inputs/logs, fit summaries and compressed importance weights are retained
so a failed stage can restart without simulating a new genome. Pathfinder CSVs
are removed after their summaries are safely saved: retaining them for 4,320
fits can require terabytes. Set `keep_stan_csv: true` if full draws are needed.
Failed-fit CSVs and MCMC chains are retained. Run folders are ignored by Git. Specify
`--config my_config.json --out runs/my_experiment` for another design; a
manifest guard prevents mixing code/config versions in a single run.

`--scenarios`, `--sources`, and `--variants` select subsets. `--only-graphs`
and `--block-limit` support debugging; partial graph sets are never reported
as complete topology comparisons. Finish all 50 block summaries before fitting.

## Outputs

`runs/default/manifest.json` records the full configuration, code/model hashes,
candidate graph/event orders, and every selected chromosome list.
`<scenario>/pool_NNN/fits/rep_NNN/<source>/<variant>/<candidate>/result.json`
contains scores, diagnostics, parameter intervals and fitted-data residuals.

After `summarize`, `summary/` contains:

- `report.md` and `topology_recovery.png`;
- `topology_recovery.csv`, `topology_rankings.csv`, `paired_ne_comparison.csv`;
- `parameter_recovery.json`, `parameter_metrics.csv`, `fit_diagnostics.csv`.

The raw tables retain pool/subsample IDs. No independent-simulation uncertainty
is inferred from reused blocks, and no best model is asserted before running
the study.

## Verification

```sh
python -B -m unittest discover -s . -p 'test_study.py' -v
```

Unit tests check topology enumeration/event ordering, haploid/diploid sizes,
sampling/exposure, strict-MRCA segmentation including path changes, raw count
aggregation, SNP ratio and jackknife SE, all candidate Stan-data shapes, HBD
inclusion, and weighted parameter summaries.
