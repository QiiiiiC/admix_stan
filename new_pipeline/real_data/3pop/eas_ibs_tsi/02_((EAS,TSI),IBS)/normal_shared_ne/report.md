# `((EAS,TSI),IBS)`

**Normal, shared Ne** | topology 02 of 21 | tree (no admixture) | 2 events, 5 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -41153.00 | +- 0.11 (MC) |
| logZ (importance sampling) | -41149.00 | |
| ESS of the IS weights | 3.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -12.13 | already applied |
| mode kept / MAP start / runtime | 2 / 2 | 5 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + TSI -> n1 | 161.6 +- 0.8 | 161.6 |
| 2 | MERGE | n1 + IBS -> root | 4.4 +- 0.5 | 166.0 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,610,161 | 0.02 |
| `IBS` | 219,437 | 0.02 |
| `TSI` | 352,882 | 0.03 |
| `n1` | 9,710 | 0.13 |
| `root` | 28,871 | 0.03 |

log-Ne random-walk step scale tau = 1.485

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,492.8 | 222 | 51.77 |
| SNP | -39,585.5 | 6 | 13209.71 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +117.33 | -115.12 | -116.79 |
| **IBS** | -115.12 | +110.65 | +116.39 |
| **TSI** | -116.79 | +116.39 | +113.19 |

![spectrum](spectrum_fit.png)
