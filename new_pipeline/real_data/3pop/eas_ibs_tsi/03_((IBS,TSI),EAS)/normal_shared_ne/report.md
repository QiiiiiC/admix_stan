# `((IBS,TSI),EAS)`

**Normal, shared Ne** | topology 03 of 21 | tree (no admixture) | 2 events, 5 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1563.14 | +- 0.22 (MC) |
| logZ (importance sampling) | -1556.45 | |
| ESS of the IS weights | 20.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -12.13 | already applied |
| mode kept / MAP start / runtime | 1 / 3 | 10 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | IBS + TSI -> n1 | 212.4 +- 0.4 | 212.4 |
| 2 | MERGE | EAS + n1 -> root | 193.5 +- 1.6 | 405.9 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 764,716 | 0.01 |
| `IBS` | 375,687 | 0.03 |
| `TSI` | 480,251 | 0.02 |
| `n1` | 1,082 | 0.01 |
| `root` | 1 | 0.24 |

log-Ne random-walk step scale tau = 1.700

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,470.2 | 222 | 52.83 |
| SNP | +28.7 | 6 | 4.99 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.94 | -0.81 | -1.06 |
| **IBS** | -0.81 | +2.99 | -1.56 |
| **TSI** | -1.06 | -1.56 | +3.48 |

![spectrum](spectrum_fit.png)
