# `((EAS,IBS),TSI)`

**Normal, separate IBD/SNP Ne** | topology 01 of 21 | tree (no admixture) | 2 events, 5 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1566.05 | +- 0.84 (MC) |
| logZ (importance sampling) | -1515.10 | |
| ESS of the IS weights | 4.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -15.06 | already applied |
| seed kept / runtime | 7 | 10 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + IBS -> n1 | 156.3 +- 1.1 | 156.3 |
| 2 | MERGE | n1 + TSI -> root | 3.6 +- 0.1 | 159.9 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,595,668 | 0.07 |
| `IBS` | 244,747 | 0.07 |
| `TSI` | 347,393 | 0.05 |
| `n1` | 10,618 | 0.02 |
| `root` | 34,788 | 0.03 |

log-Ne random-walk step scale tau_ibd = 1.559

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 867 | 0.02 |
| `IBS` | 121,586 | 0.08 |
| `TSI` | 128,767 | 0.07 |
| `n1` | 11,567 | 0.02 |
| `root` | 14,195 | 0.02 |

log-Ne random-walk step scale tau_snp = 0.884

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,460.1 | 222 | 51.39 |
| SNP | +32.9 | 6 | 3.58 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.09 | -0.15 | -0.02 |
| **IBS** | -0.15 | +0.02 | +0.30 |
| **TSI** | -0.02 | +0.30 | -0.23 |

![spectrum](spectrum_fit.png)
