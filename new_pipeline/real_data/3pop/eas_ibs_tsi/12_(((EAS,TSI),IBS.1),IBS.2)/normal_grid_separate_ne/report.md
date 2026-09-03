# `(((EAS,TSI),IBS.1),IBS.2)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 12 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1062.75 | +- 0.34 (MC) |
| logZ (importance sampling) | -1033.65 | |
| ESS of the IS weights | 1.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 1 | 25 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 85.2 +- 1.0 | 85.2 |
| 2 | MERGE | EAS + TSI -> n1 | 71.5 +- 0.4 | 156.7 |
| 3 | MERGE | IBS.2 + n1 -> n2 | 2.3 +- 0.0 | 159.0 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.4 +- 0.0 | 160.4 |

## Admixture fraction

**f = 0.994 +- 0.000** (fraction from `IBS.1`; 0.006 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 25,059,675 | 18,105,172 |
| `IBS` | 2,134,275 | 1,732,810 |
| `TSI` | 1,192,553 | 946,315 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,355,339 | 0.03 |
| `IBS` | 743,630 | 0.09 |
| `TSI` | 347,059 | 0.02 |
| `IBS.1` | 46,140 | 0.02 |
| `IBS.2` | 21,055 | 0.01 |
| `n1` | 16,085 | 0.01 |
| `n2` | 21,122 | 0.01 |
| `root` | 34,194 | 0.01 |

log-Ne random-walk step scale tau_ibd = 2.076

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 815 | 829 |
| `IBS` | 196,938 | 197,945 |
| `TSI` | 108,548 | 110,382 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 872 | 0.01 |
| `IBS` | 204,773 | 0.05 |
| `TSI` | 114,834 | 0.03 |
| `IBS.1` | 81,070 | 0.05 |
| `IBS.2` | 31,662 | 0.01 |
| `n1` | 29,772 | 0.01 |
| `n2` | 31,617 | 0.01 |
| `root` | 32,975 | 0.01 |

log-Ne random-walk step scale tau_snp = 1.130

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -821.2 | 222 | 45.98 |
| SNP | +35.0 | 6 | 2.87 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.03 | +0.13 | -0.08 |
| **IBS** | +0.13 | -0.39 | +0.14 |
| **TSI** | -0.08 | +0.14 | +0.01 |

![spectrum](spectrum_fit.png)
