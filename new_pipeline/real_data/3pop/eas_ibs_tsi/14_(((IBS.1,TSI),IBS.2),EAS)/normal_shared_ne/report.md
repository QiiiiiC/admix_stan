# `(((IBS.1,TSI),IBS.2),EAS)`

**Normal, shared Ne** | topology 14 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +31.52 | +- 0.34 (MC) |
| logZ (importance sampling) | +53.93 | |
| ESS of the IS weights | 12.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 1 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 64.3 +- 0.8 | 64.3 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 61.7 +- 0.2 | 126.0 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 202.9 +- 0.8 | 328.9 |
| 4 | MERGE | n2 + EAS -> root | 16.1 +- 0.0 | 344.9 |

## Admixture fraction

**f = 0.094 +- 0.001** (fraction from `IBS.1`; 0.906 from `IBS.2`)

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 964,120 | 0.02 |
| `IBS` | 1,014,687 | 0.07 |
| `TSI` | 402,782 | 0.07 |
| `IBS.1` | 1,288 | 0.02 |
| `IBS.2` | 184,931 | 0.03 |
| `n1` | 69,576 | 0.03 |
| `n2` | 90 | 0.01 |
| `root` | 50 | 0.01 |

log-Ne random-walk step scale tau = 1.443

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +155.3 | 222 | 37.71 |
| SNP | +34.4 | 6 | 3.07 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.09 | +0.78 | -0.61 |
| **IBS** | +0.78 | -1.10 | -0.42 |
| **TSI** | -0.61 | -0.42 | +1.55 |

![spectrum](spectrum_fit.png)
