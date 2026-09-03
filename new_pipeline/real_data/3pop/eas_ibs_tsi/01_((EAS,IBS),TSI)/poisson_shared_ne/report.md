# `((EAS,IBS),TSI)`

**Poisson, shared Ne** | topology 01 of 21 | tree (no admixture) | 2 events, 5 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -43568.24 | +- 0.29 (MC) |
| logZ (importance sampling) | -43559.95 | |
| ESS of the IS weights | 2.9 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -12.13 | already applied |
| seed kept / runtime | 7 | 7 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + IBS -> n1 | 283.1 +- 0.8 | 283.1 |
| 2 | MERGE | n1 + TSI -> root | 7.8 +- 0.4 | 290.9 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,185,340 | 0.02 |
| `IBS` | 241,230 | 0.02 |
| `TSI` | 326,021 | 0.02 |
| `n1` | 271 | 0.05 |
| `root` | 721 | 0.03 |

log-Ne random-walk step scale tau = 1.545

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,576.9 | 222 | 6327.28 |
| SNP | -35,897.0 | 6 | 11980.21 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +112.52 | -124.35 | -97.96 |
| **IBS** | -124.35 | +92.95 | +154.21 |
| **TSI** | -97.96 | +154.21 | +41.70 |

![spectrum](spectrum_fit.png)
