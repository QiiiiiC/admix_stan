# `((EAS,IBS),TSI)`

**Poisson, separate IBD/SNP Ne** | topology 01 of 21 | tree (no admixture) | 2 events, 5 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -7618.19 | +- 0.20 (MC) |
| logZ (importance sampling) | -7602.00 | |
| ESS of the IS weights | 4.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -15.06 | already applied |
| mode kept / MAP start / runtime | 1 / 11 | 11 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + IBS -> n1 | 228.9 +- 0.8 | 228.9 |
| 2 | MERGE | n1 + TSI -> root | 11.3 +- 0.3 | 240.2 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,225,218 | 0.03 |
| `IBS` | 253,291 | 0.02 |
| `TSI` | 329,656 | 0.03 |
| `n1` | 2,702 | 0.03 |
| `root` | 4,653 | 0.03 |

log-Ne random-walk step scale tau_ibd = 1.260

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,269 | 0.01 |
| `IBS` | 146,150 | 0.10 |
| `TSI` | 329,213 | 0.06 |
| `n1` | 21,722 | 0.03 |
| `root` | 26,770 | 0.03 |

log-Ne random-walk step scale tau_snp = 0.862

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,525.9 | 222 | 384.26 |
| SNP | +39.4 | 6 | 1.39 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.00 | +0.22 | -0.22 |
| **IBS** | +0.22 | -0.58 | +0.18 |
| **TSI** | -0.22 | +0.18 | +0.26 |

![spectrum](spectrum_fit.png)
