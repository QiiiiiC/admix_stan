# `((IBS,TSI),EAS)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 03 of 21 | tree (no admixture) | 2 events, 5 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +54.47 | +- 1.15 (MC) |
| logZ (importance sampling) | +118.96 | |
| ESS of the IS weights | 1.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -26.08 | already applied |
| seed kept / runtime | 13 | 15 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | IBS + TSI -> n1 | 147.4 +- 1.7 | 147.4 |
| 2 | MERGE | EAS + n1 -> root | 192.9 +- 8.5 | 340.3 |

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 26,128,272 | 16,343,564 |
| `IBS` | 2,458,237 | 1,707,613 |
| `TSI` | 1,268,389 | 886,569 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 865,326 | 0.03 |
| `IBS` | 278,776 | 0.10 |
| `TSI` | 404,181 | 0.06 |
| `n1` | 24,816 | 0.08 |
| `root` | 93 | 0.38 |

log-Ne random-walk step scale tau_ibd = 1.978

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 2,269 | 2,211 |
| `IBS` | 107,377 | 104,590 |
| `TSI` | 144,522 | 142,524 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,928 | 0.03 |
| `IBS` | 93,518 | 0.04 |
| `TSI` | 142,326 | 0.04 |
| `n1` | 41,153 | 0.05 |
| `root` | 12,550 | 0.08 |

log-Ne random-walk step scale tau_snp = 0.405

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +293.6 | 222 | 36.64 |
| SNP | +28.8 | 6 | 4.94 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.08 | +0.39 | -0.24 |
| **IBS** | +0.39 | -0.52 | -0.24 |
| **TSI** | -0.24 | -0.24 | +0.69 |

![spectrum](spectrum_fit.png)
