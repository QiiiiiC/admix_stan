# `(((EAS,TSI.1),TSI.2),IBS)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 18 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +89.64 | +- 0.35 (MC) |
| logZ (importance sampling) | +112.47 | |
| ESS of the IS weights | 10.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 7 | 25 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 113.1 +- 0.4 | 113.1 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 17.8 +- 0.2 | 130.9 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 15.6 +- 0.3 | 146.5 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 147.6 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,291,201 | 3,651,895 |
| `IBS` | 2,397,862 | 1,664,509 |
| `TSI` | 1,355,973 | 964,029 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,546,644 | 0.06 |
| `IBS` | 211,559 | 0.03 |
| `TSI` | 393,908 | 0.06 |
| `TSI.1` | 61,956 | 0.03 |
| `TSI.2` | 21,946 | 0.02 |
| `n1` | 23,743 | 0.02 |
| `n2` | 158,481 | 0.02 |
| `root` | 74,872 | 0.01 |

log-Ne random-walk step scale tau_ibd = 2.028

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 764 | 765 |
| `IBS` | 91,232 | 94,728 |
| `TSI` | 125,731 | 124,720 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 730 | 0.01 |
| `IBS` | 96,840 | 0.06 |
| `TSI` | 123,889 | 0.06 |
| `TSI.1` | 41,756 | 0.05 |
| `TSI.2` | 11,258 | 0.02 |
| `n1` | 11,398 | 0.02 |
| `n2` | 17,509 | 0.03 |
| `root` | 17,450 | 0.03 |

log-Ne random-walk step scale tau_snp = 1.025

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +337.5 | 222 | 35.34 |
| SNP | +36.7 | 6 | 2.29 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.08 | -0.07 | -0.08 |
| **IBS** | -0.07 | -0.75 | +0.95 |
| **TSI** | -0.08 | +0.95 | -0.75 |

![spectrum](spectrum_fit.png)
