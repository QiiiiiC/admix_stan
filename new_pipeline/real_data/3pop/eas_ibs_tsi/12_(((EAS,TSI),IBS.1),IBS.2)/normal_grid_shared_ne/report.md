# `(((EAS,TSI),IBS.1),IBS.2)`

**Normal, shared Ne, recent grid** | topology 12 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -37200.12 | +- 0.35 (MC) |
| logZ (importance sampling) | -37174.38 | |
| ESS of the IS weights | 2.8 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 1 / 7 | 30 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 110.2 +- 0.5 | 110.2 |
| 2 | MERGE | EAS + TSI -> n1 | 46.2 +- 1.3 | 156.4 |
| 3 | MERGE | IBS.2 + n1 -> n2 | 4.6 +- 0.3 | 161.0 |
| 4 | MERGE | IBS.1 + n2 -> root | 2.1 +- 0.1 | 163.1 |

## Admixture fraction

**f = 0.035 +- 0.001** (fraction from `IBS.1`; 0.965 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 12,736,106 | 13,467,536 |
| `IBS` | 1,711,637 | 1,313,617 |
| `TSI` | 855,036 | 701,352 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,316,347 | 0.03 |
| `IBS` | 715,561 | 0.12 |
| `TSI` | 339,106 | 0.05 |
| `IBS.1` | 2 | 0.07 |
| `IBS.2` | 161,602 | 0.04 |
| `n1` | 30,067 | 0.06 |
| `n2` | 47,865 | 0.05 |
| `root` | 32,061 | 0.02 |

log-Ne random-walk step scale tau = 1.759

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -841.4 | 222 | 46.15 |
| SNP | -36,159.2 | 6 | 12067.62 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +112.74 | -97.06 | -125.87 |
| **IBS** | -97.06 | +41.91 | +153.20 |
| **TSI** | -125.87 | +153.20 | +95.88 |

![spectrum](spectrum_fit.png)
