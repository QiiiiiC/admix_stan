# `((IBS,TSI),EAS)`

**Poisson, separate IBD/SNP Ne** | topology 03 of 21 | tree (no admixture) | 2 events, 5 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5783.41 | +- 0.12 (MC) |
| logZ (importance sampling) | -5774.53 | |
| ESS of the IS weights | 18.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -15.06 | already applied |
| seed kept / runtime | 13 | 9 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | IBS + TSI -> n1 | 160.7 +- 1.0 | 160.7 |
| 2 | MERGE | EAS + n1 -> root | 151.1 +- 6.5 | 311.8 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 895,963 | 0.02 |
| `IBS` | 314,650 | 0.04 |
| `TSI` | 424,309 | 0.03 |
| `n1` | 13,854 | 0.03 |
| `root` | 887 | 0.29 |

log-Ne random-walk step scale tau_ibd = 1.332

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 209,292 | 0.26 |
| `IBS` | 126,735 | 0.09 |
| `TSI` | 101,617 | 0.12 |
| `n1` | 846 | 0.04 |
| `root` | 11,647 | 0.01 |

log-Ne random-walk step scale tau_snp = 1.312

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,673.8 | 222 | 46.40 |
| SNP | +41.8 | 6 | 0.59 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.05 | -0.15 | +0.04 |
| **IBS** | -0.15 | -0.01 | +0.30 |
| **TSI** | +0.04 | +0.30 | -0.36 |

![spectrum](spectrum_fit.png)
