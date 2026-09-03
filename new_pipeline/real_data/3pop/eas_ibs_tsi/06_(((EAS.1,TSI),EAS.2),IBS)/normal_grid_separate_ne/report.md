# `(((EAS.1,TSI),EAS.2),IBS)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 06 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +53.54 | +- 0.17 (MC) |
| logZ (importance sampling) | +71.47 | |
| ESS of the IS weights | 1.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 1 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 131.5 +- 0.4 | 131.5 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 13.9 +- 0.1 | 145.4 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 1.0 +- 0.0 | 146.4 |
| 4 | MERGE | n2 + IBS -> root | 1.1 +- 0.0 | 147.5 |

## Admixture fraction

**f = 0.012 +- 0.000** (fraction from `EAS.1`; 0.988 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 3,759,776 | 3,365,294 |
| `IBS` | 3,535,366 | 2,294,620 |
| `TSI` | 2,521,278 | 1,680,373 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,572,637 | 0.02 |
| `IBS` | 211,177 | 0.02 |
| `TSI` | 317,101 | 0.04 |
| `EAS.1` | 267,787 | 0.03 |
| `EAS.2` | 20,653 | 0.01 |
| `n1` | 354,663 | 0.02 |
| `n2` | 250,046 | 0.01 |
| `root` | 74,834 | 0.01 |

log-Ne random-walk step scale tau_ibd = 2.497

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 730 | 732 |
| `IBS` | 118,016 | 116,978 |
| `TSI` | 101,923 | 101,251 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 734 | 0.01 |
| `IBS` | 112,159 | 0.08 |
| `TSI` | 98,418 | 0.07 |
| `EAS.1` | 24,629 | 0.02 |
| `EAS.2` | 13,294 | 0.02 |
| `n1` | 23,428 | 0.01 |
| `n2` | 23,921 | 0.01 |
| `root` | 24,066 | 0.01 |

log-Ne random-walk step scale tau_snp = 1.168

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +317.9 | 222 | 35.50 |
| SNP | +39.1 | 6 | 1.51 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.05 | -0.02 | -0.08 |
| **IBS** | -0.02 | -0.24 | +0.31 |
| **TSI** | -0.08 | +0.31 | -0.14 |

![spectrum](spectrum_fit.png)
