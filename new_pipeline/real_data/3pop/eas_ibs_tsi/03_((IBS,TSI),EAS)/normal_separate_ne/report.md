# `((IBS,TSI),EAS)`

**Normal, separate IBD/SNP Ne** | topology 03 of 21 | tree (no admixture) | 2 events, 5 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -334.32 | +- 1.12 (MC) |
| logZ (importance sampling) | -266.34 | |
| ESS of the IS weights | 2.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -15.06 | already applied |
| seed kept / runtime | 1 | 10 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | IBS + TSI -> n1 | 143.2 +- 0.8 | 143.2 |
| 2 | MERGE | EAS + n1 -> root | 205.4 +- 9.3 | 348.6 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 967,768 | 0.03 |
| `IBS` | 286,158 | 0.08 |
| `TSI` | 419,925 | 0.06 |
| `n1` | 30,158 | 0.04 |
| `root` | 43 | 0.53 |

log-Ne random-walk step scale tau_ibd = 1.526

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 72,623 | 0.18 |
| `IBS` | 107,994 | 0.05 |
| `TSI` | 113,656 | 0.19 |
| `n1` | 1,169 | 0.05 |
| `root` | 5,541 | 0.19 |

log-Ne random-walk step scale tau_snp = 1.233

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -144.5 | 222 | 40.37 |
| SNP | +7.6 | 6 | 12.00 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.22 | +0.29 | +0.14 |
| **IBS** | +0.29 | -0.24 | -0.33 |
| **TSI** | +0.14 | -0.33 | +0.04 |

![spectrum](spectrum_fit.png)
