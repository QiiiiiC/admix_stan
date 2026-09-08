# `((EAS,TSI),IBS)`

**Poisson, separate IBD/SNP Ne** | topology 02 of 21 | tree (no admixture) | 2 events, 5 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -7683.70 | +- 0.26 (MC) |
| logZ (importance sampling) | -7665.43 | |
| ESS of the IS weights | 1.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -15.06 | already applied |
| mode kept / MAP start / runtime | 1 / 9 | 11 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + TSI -> n1 | 229.1 +- 0.5 | 229.1 |
| 2 | MERGE | n1 + IBS -> root | 8.7 +- 0.4 | 237.8 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,195,580 | 0.02 |
| `IBS` | 243,076 | 0.04 |
| `TSI` | 343,900 | 0.03 |
| `n1` | 2,679 | 0.02 |
| `root` | 4,808 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.349

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,271 | 0.01 |
| `IBS` | 221,857 | 0.07 |
| `TSI` | 166,774 | 0.06 |
| `n1` | 20,022 | 0.01 |
| `root` | 22,711 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.930

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,592.8 | 222 | 328.23 |
| SNP | +37.5 | 6 | 2.05 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.02 | +0.10 | -0.15 |
| **IBS** | +0.10 | -0.46 | +0.28 |
| **TSI** | -0.15 | +0.28 | +0.02 |

![spectrum](spectrum_fit.png)
