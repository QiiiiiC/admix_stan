# `((EAS,IBS.1),(IBS.2,TSI))`

**Normal, shared Ne** | topology 15 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +3079.62 | +- 0.49 (MC) |
| logZ (importance sampling) | +3108.46 | |
| ESS of the IS weights | 1.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| mode kept / MAP start / runtime | 2 / 2 | 21 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 51.8 +- 0.3 | 51.8 |
| 2 | MERGE | IBS.1 + EAS -> n1 | 149.4 +- 0.3 | 201.2 |
| 3 | MERGE | IBS.2 + TSI -> n2 | 1.1 +- 0.0 | 202.4 |
| 4 | MERGE | n1 + n2 -> root | 246.1 +- 1.0 | 448.5 |

## Admixture fraction

**f = 0.006 +- 0.000** (fraction from `IBS.1`; 0.994 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,956,420 | 0.04 |
| `IBS` | 1,140,449 | 0.06 |
| `TSI` | 453,539 | 0.07 |
| `IBS.1` | 7,288 | 0.04 |
| `IBS.2` | 146,598 | 0.03 |
| `n1` | 3,981 | 0.03 |
| `n2` | 2,092 | 0.01 |
| `root` | 13 | 0.00 |

log-Ne random-walk step scale tau = 1.302

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +3,281.7 | 222 | 10.20 |
| SNP | +26.4 | 6 | 5.75 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +1.24 | -1.94 | -0.51 |
| **IBS** | -1.94 | +4.16 | -0.51 |
| **TSI** | -0.51 | -0.51 | +1.46 |

![spectrum](spectrum_fit.png)
