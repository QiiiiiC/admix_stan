# `(((EAS,TSI.1),TSI.2),IBS)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 18 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -4324.86 | +- 0.98 (MC) |
| logZ (importance sampling) | -4254.99 | |
| ESS of the IS weights | 1.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 13 | 22 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 140.9 +- 1.7 | 140.9 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 1.0 +- 0.0 | 141.9 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 34.9 +- 0.7 | 176.8 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 177.8 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 3,771,073 | 3,632,700 |
| `IBS` | 3,785,127 | 2,367,966 |
| `TSI` | 2,297,052 | 1,543,478 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,387,737 | 0.06 |
| `IBS` | 223,349 | 0.04 |
| `TSI` | 304,551 | 0.09 |
| `TSI.1` | 535,193 | 0.09 |
| `TSI.2` | 22,463 | 0.03 |
| `n1` | 22,388 | 0.03 |
| `n2` | 146,685 | 0.02 |
| `root` | 75,461 | 0.01 |

log-Ne random-walk step scale tau_ibd = 2.106

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,326 | 1,217 |
| `IBS` | 258,110 | 237,177 |
| `TSI` | 493,046 | 468,963 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 878 | 0.02 |
| `IBS` | 170,287 | 0.13 |
| `TSI` | 393,531 | 0.28 |
| `TSI.1` | 23,599 | 0.10 |
| `TSI.2` | 1,585 | 0.02 |
| `n1` | 1,583 | 0.02 |
| `n2` | 9,862 | 0.03 |
| `root` | 8,690 | 0.03 |

log-Ne random-walk step scale tau_snp = 0.896

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,048.1 | 222 | 77.45 |
| SNP | +24.6 | 6 | 6.35 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.06 | -0.37 | +0.25 |
| **IBS** | -0.37 | +0.16 | +0.57 |
| **TSI** | +0.25 | +0.57 | -1.01 |

![spectrum](spectrum_fit.png)
