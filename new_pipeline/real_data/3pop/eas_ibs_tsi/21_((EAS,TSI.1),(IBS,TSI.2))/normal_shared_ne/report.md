# `((EAS,TSI.1),(IBS,TSI.2))`

**Normal, shared Ne** | topology 21 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -42829.08 | +- 3.69 (MC) |
| logZ (importance sampling) | -42605.76 | |
| ESS of the IS weights | 1.9 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 13 | 10 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | TSI.1 + EAS -> n1 | 1.0 +- 0.0 | 2.0 |
| 3 | MERGE | TSI.2 + IBS -> n2 | 1.0 +- 0.0 | 3.0 |
| 4 | MERGE | n1 + n2 -> root | 340.2 +- 8.9 | 343.2 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `TSI.1`; 1.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,091,203 | 0.13 |
| `IBS` | 366,265 | 0.16 |
| `TSI` | 361,589 | 0.16 |
| `TSI.1` | 1,002,287 | 0.08 |
| `TSI.2` | 362,445 | 0.16 |
| `n1` | 1,002,307 | 0.08 |
| `n2` | 360,982 | 0.16 |
| `root` | 46 | 0.79 |

log-Ne random-walk step scale tau = 1.520

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -3,104.5 | 222 | 64.02 |
| SNP | -39,467.2 | 6 | 13170.28 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +116.81 | -114.97 | -115.91 |
| **IBS** | -114.97 | +113.12 | +113.44 |
| **TSI** | -115.91 | +113.44 | +114.29 |

![spectrum](spectrum_fit.png)
