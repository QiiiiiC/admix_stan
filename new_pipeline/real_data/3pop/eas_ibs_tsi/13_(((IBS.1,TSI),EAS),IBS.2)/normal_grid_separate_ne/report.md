# `(((IBS.1,TSI),EAS),IBS.2)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 13 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +779.51 | +- 0.62 (MC) |
| logZ (importance sampling) | +820.81 | |
| ESS of the IS weights | 1.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 1 | 24 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 87.3 +- 2.6 | 87.3 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 31.9 +- 0.9 | 119.2 |
| 3 | MERGE | n1 + EAS -> n2 | 237.2 +- 2.6 | 356.4 |
| 4 | MERGE | IBS.1 + n2 -> root | 365.0 +- 18.4 | 721.4 |

## Admixture fraction

**f = 0.689 +- 0.007** (fraction from `IBS.1`; 0.311 from `IBS.2`)

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 26,647,172 | 26,162,994 |
| `IBS` | 665,974 | 725,072 |
| `TSI` | 3,863,758 | 2,537,398 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,075,579 | 0.05 |
| `IBS` | 884,571 | 0.08 |
| `TSI` | 495,959 | 0.08 |
| `IBS.1` | 44,875 | 0.07 |
| `IBS.2` | 5,026 | 0.10 |
| `n1` | 29,717 | 0.08 |
| `n2` | 6 | 0.51 |
| `root` | 12,138 | 0.00 |

log-Ne random-walk step scale tau_ibd = 2.515

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,474 | 1,586 |
| `IBS` | 550,780 | 567,435 |
| `TSI` | 1,676,214 | 1,696,565 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,993 | 0.02 |
| `IBS` | 638,346 | 0.17 |
| `TSI` | 1,696,730 | 0.18 |
| `IBS.1` | 189,494 | 0.17 |
| `IBS.2` | 8,914,131 | 0.25 |
| `n1` | 2,953,370 | 0.20 |
| `n2` | 330,611 | 0.15 |
| `root` | 109,908 | 0.12 |

log-Ne random-walk step scale tau_snp = 0.313

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +1,118.1 | 222 | 29.16 |
| SNP | +27.6 | 6 | 5.33 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.02 | +1.09 | -1.13 |
| **IBS** | +1.09 | -1.69 | -0.40 |
| **TSI** | -1.13 | -0.40 | +2.53 |

![spectrum](spectrum_fit.png)
