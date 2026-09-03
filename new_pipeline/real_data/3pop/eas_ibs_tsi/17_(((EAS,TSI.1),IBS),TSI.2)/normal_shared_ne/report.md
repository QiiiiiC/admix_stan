# `(((EAS,TSI.1),IBS),TSI.2)`

**Normal, shared Ne** | topology 17 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -49362.66 | +- 0.39 (MC) |
| logZ (importance sampling) | -49339.06 | |
| ESS of the IS weights | 2.8 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 7 | 9 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 1.0 +- 0.0 | 2.0 |
| 3 | MERGE | n1 + IBS -> n2 | 1.0 +- 0.0 | 3.0 |
| 4 | MERGE | TSI.1 + n2 -> root | 150.8 +- 1.0 | 153.8 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,877,623 | 0.06 |
| `IBS` | 1,445,870 | 0.06 |
| `TSI` | 358,084 | 0.04 |
| `TSI.1` | 359,582 | 0.04 |
| `TSI.2` | 1,761,026 | 0.06 |
| `n1` | 1,761,063 | 0.06 |
| `n2` | 1,596,219 | 0.06 |
| `root` | 25,908 | 0.02 |

log-Ne random-walk step scale tau = 0.918

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -9,257.7 | 222 | 119.06 |
| SNP | -39,911.6 | 6 | 13318.41 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +117.58 | -115.98 | -116.42 |
| **IBS** | -115.98 | +113.62 | +114.96 |
| **TSI** | -116.42 | +114.96 | +113.83 |

![spectrum](spectrum_fit.png)
