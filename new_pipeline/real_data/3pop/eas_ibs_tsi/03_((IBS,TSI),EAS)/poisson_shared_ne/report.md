# `((IBS,TSI),EAS)`

**Poisson, shared Ne** | topology 03 of 21 | tree (no admixture) | 2 events, 5 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5834.21 | +- 0.30 (MC) |
| logZ (importance sampling) | -5825.54 | |
| ESS of the IS weights | 13.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -12.13 | already applied |
| mode kept / MAP start / runtime | 1 / 3 | 10 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | IBS + TSI -> n1 | 238.8 +- 0.8 | 238.8 |
| 2 | MERGE | EAS + n1 -> root | 151.2 +- 2.5 | 390.0 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 892,530 | 0.02 |
| `IBS` | 312,045 | 0.03 |
| `TSI` | 418,345 | 0.03 |
| `n1` | 837 | 0.02 |
| `root` | 21 | 0.18 |

log-Ne random-walk step scale tau = 1.490

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,768.5 | 222 | 96.84 |
| SNP | +34.6 | 6 | 3.01 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.20 | +0.35 | +0.04 |
| **IBS** | +0.35 | +1.40 | -2.22 |
| **TSI** | +0.04 | -2.22 | +2.01 |

![spectrum](spectrum_fit.png)
