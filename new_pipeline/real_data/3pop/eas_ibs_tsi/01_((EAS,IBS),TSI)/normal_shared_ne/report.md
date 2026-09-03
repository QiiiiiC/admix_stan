# `((EAS,IBS),TSI)`

**Normal, shared Ne** | topology 01 of 21 | tree (no admixture) | 2 events, 5 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -41110.24 | +- 0.03 (MC) |
| logZ (importance sampling) | -41106.95 | |
| ESS of the IS weights | 4.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -12.13 | already applied |
| seed kept / runtime | 7 | 5 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + IBS -> n1 | 161.8 +- 0.5 | 161.8 |
| 2 | MERGE | n1 + TSI -> root | 4.2 +- 0.5 | 166.0 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,645,380 | 0.02 |
| `IBS` | 235,568 | 0.02 |
| `TSI` | 328,477 | 0.02 |
| `n1` | 7,922 | 0.12 |
| `root` | 29,632 | 0.02 |

log-Ne random-walk step scale tau = 1.513

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,461.8 | 222 | 51.49 |
| SNP | -39,568.7 | 6 | 13204.11 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +117.32 | -115.64 | -116.24 |
| **IBS** | -115.64 | +111.57 | +116.47 |
| **TSI** | -116.24 | +116.47 | +112.08 |

![spectrum](spectrum_fit.png)
