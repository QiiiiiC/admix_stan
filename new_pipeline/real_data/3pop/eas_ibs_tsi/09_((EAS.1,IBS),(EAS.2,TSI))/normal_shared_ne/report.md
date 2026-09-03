# `((EAS.1,IBS),(EAS.2,TSI))`

**Normal, shared Ne** | topology 09 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -39242.55 | +- 0.87 (MC) |
| logZ (importance sampling) | -39180.16 | |
| ESS of the IS weights | 7.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 13 | 19 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 171.0 +- 2.3 | 171.0 |
| 2 | MERGE | EAS.1 + IBS -> n1 | 3.7 +- 1.0 | 174.7 |
| 3 | MERGE | EAS.2 + TSI -> n2 | 1.0 +- 0.0 | 175.7 |
| 4 | MERGE | n1 + n2 -> root | 1.0 +- 0.0 | 176.7 |

## Admixture fraction

**f = 0.999 +- 0.000** (fraction from `EAS.1`; 0.001 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,492,254 | 0.06 |
| `IBS` | 213,399 | 0.09 |
| `TSI` | 320,427 | 0.09 |
| `EAS.1` | 1,664 | 0.35 |
| `EAS.2` | 14,488 | 0.31 |
| `n1` | 143,486 | 0.07 |
| `n2` | 14,479 | 0.32 |
| `root` | 27,008 | 0.04 |

log-Ne random-walk step scale tau = 3.317

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -352.8 | 222 | 41.84 |
| SNP | -38,651.9 | 6 | 12898.52 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +115.90 | -113.98 | -115.09 |
| **IBS** | -113.98 | +110.11 | +114.66 |
| **TSI** | -115.09 | +114.66 | +111.59 |

![spectrum](spectrum_fit.png)
