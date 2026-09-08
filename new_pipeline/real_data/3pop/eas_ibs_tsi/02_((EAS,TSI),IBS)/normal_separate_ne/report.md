# `((EAS,TSI),IBS)`

**Normal, separate IBD/SNP Ne** | topology 02 of 21 | tree (no admixture) | 2 events, 5 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1545.83 | +- 0.12 (MC) |
| logZ (importance sampling) | -1536.81 | |
| ESS of the IS weights | 4.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -15.06 | already applied |
| mode kept / MAP start / runtime | 1 / 8 | 11 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + TSI -> n1 | 155.6 +- 0.3 | 155.6 |
| 2 | MERGE | n1 + IBS -> root | 5.1 +- 0.1 | 160.7 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,540,563 | 0.02 |
| `IBS` | 231,636 | 0.03 |
| `TSI` | 365,319 | 0.03 |
| `n1` | 18,602 | 0.01 |
| `root` | 32,737 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.266

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 863 | 0.01 |
| `IBS` | 145,250 | 0.09 |
| `TSI` | 106,631 | 0.08 |
| `n1` | 17,622 | 0.01 |
| `root` | 18,460 | 0.01 |

log-Ne random-walk step scale tau_snp = 1.006

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,465.2 | 222 | 51.44 |
| SNP | +41.5 | 6 | 0.72 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.03 | +0.00 | -0.06 |
| **IBS** | +0.00 | -0.25 | +0.26 |
| **TSI** | -0.06 | +0.26 | -0.14 |

![spectrum](spectrum_fit.png)
