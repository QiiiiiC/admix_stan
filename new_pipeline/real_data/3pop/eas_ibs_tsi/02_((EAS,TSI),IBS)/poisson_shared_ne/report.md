# `((EAS,TSI),IBS)`

**Poisson, shared Ne** | topology 02 of 21 | tree (no admixture) | 2 events, 5 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -43992.11 | +- 0.04 (MC) |
| logZ (importance sampling) | -43989.58 | |
| ESS of the IS weights | 56.6 / 4000 | ok |
| Stan dropped-constant correction | -12.13 | already applied |
| seed kept / runtime | 13 | 4 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + TSI -> n1 | 292.6 +- 1.9 | 292.6 |
| 2 | MERGE | n1 + IBS -> root | 1.0 +- 0.0 | 293.6 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,169,003 | 0.01 |
| `IBS` | 241,492 | 0.02 |
| `TSI` | 328,788 | 0.02 |
| `n1` | 36 | 0.01 |
| `root` | 608 | 0.09 |

log-Ne random-walk step scale tau = 2.995

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,643.0 | 222 | 6958.32 |
| SNP | -36,043.4 | 6 | 12029.00 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +112.63 | -97.33 | -125.39 |
| **IBS** | -97.33 | +42.36 | +153.26 |
| **TSI** | -125.39 | +153.26 | +94.91 |

![spectrum](spectrum_fit.png)
