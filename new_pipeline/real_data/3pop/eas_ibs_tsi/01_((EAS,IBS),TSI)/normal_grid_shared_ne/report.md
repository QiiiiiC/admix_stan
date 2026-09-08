# `((EAS,IBS),TSI)`

**Normal, shared Ne, recent grid** | topology 01 of 21 | tree (no admixture) | 2 events, 5 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -40954.05 | +- 0.35 (MC) |
| logZ (importance sampling) | -40930.89 | |
| ESS of the IS weights | 13.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -17.65 | already applied |
| mode kept / MAP start / runtime | 1 / 6 | 10 s |
| mode search | 12/12 MAP starts succeeded | 2 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + IBS -> n1 | 163.0 +- 0.6 | 163.0 |
| 2 | MERGE | n1 + TSI -> root | 3.1 +- 0.1 | 166.1 |

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 25,465,255 | 17,673,041 |
| `IBS` | 3,609,258 | 2,237,549 |
| `TSI` | 1,933,785 | 1,423,469 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,454,736 | 0.04 |
| `IBS` | 224,233 | 0.05 |
| `TSI` | 316,883 | 0.07 |
| `n1` | 7,138 | 0.06 |
| `root` | 29,391 | 0.02 |

log-Ne random-walk step scale tau = 2.429

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,218.1 | 222 | 49.52 |
| SNP | -39,590.8 | 6 | 13211.48 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +117.33 | -115.60 | -116.31 |
| **IBS** | -115.60 | +111.63 | +116.33 |
| **TSI** | -116.31 | +116.33 | +112.33 |

![spectrum](spectrum_fit.png)
