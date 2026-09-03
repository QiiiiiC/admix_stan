# `((EAS,TSI),IBS)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 02 of 21 | tree (no admixture) | 2 events, 5 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1384.85 | +- 0.31 (MC) |
| logZ (importance sampling) | -1358.66 | |
| ESS of the IS weights | 2.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -26.08 | already applied |
| seed kept / runtime | 13 | 16 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + TSI -> n1 | 160.9 +- 0.5 | 160.9 |
| 2 | MERGE | n1 + IBS -> root | 1.0 +- 0.0 | 161.9 |

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 31,979,160 | 22,432,203 |
| `IBS` | 3,994,039 | 2,565,974 |
| `TSI` | 2,871,730 | 1,825,997 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,282,167 | 0.05 |
| `IBS` | 224,826 | 0.02 |
| `TSI` | 340,336 | 0.05 |
| `n1` | 13,663 | 0.03 |
| `root` | 30,022 | 0.03 |

log-Ne random-walk step scale tau_ibd = 2.628

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 878 | 881 |
| `IBS` | 145,262 | 145,305 |
| `TSI` | 95,700 | 95,664 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 893 | 0.01 |
| `IBS` | 145,956 | 0.07 |
| `TSI` | 95,728 | 0.06 |
| `n1` | 16,891 | 0.01 |
| `root` | 16,967 | 0.00 |

log-Ne random-walk step scale tau_snp = 1.033

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,198.1 | 222 | 49.29 |
| SNP | +38.2 | 6 | 1.79 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.00 | -0.18 | +0.18 |
| **IBS** | -0.18 | +0.15 | +0.22 |
| **TSI** | +0.18 | +0.22 | -0.55 |

![spectrum](spectrum_fit.png)
