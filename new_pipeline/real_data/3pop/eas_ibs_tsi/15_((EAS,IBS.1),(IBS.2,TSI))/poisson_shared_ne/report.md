# `((EAS,IBS.1),(IBS.2,TSI))`

**Poisson, shared Ne** | topology 15 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1013.02 | +- 0.70 (MC) |
| logZ (importance sampling) | -970.66 | |
| ESS of the IS weights | 4.9 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 1 | 15 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 58.3 +- 0.1 | 58.3 |
| 2 | MERGE | IBS.1 + EAS -> n1 | 157.4 +- 0.7 | 215.7 |
| 3 | MERGE | IBS.2 + TSI -> n2 | 1.0 +- 0.0 | 216.7 |
| 4 | MERGE | n1 + n2 -> root | 201.6 +- 0.8 | 418.3 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `IBS.1`; 1.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,894,858 | 0.04 |
| `IBS` | 1,064,893 | 0.10 |
| `TSI` | 413,077 | 0.11 |
| `IBS.1` | 2,714 | 0.02 |
| `IBS.2` | 112,262 | 0.06 |
| `n1` | 2,529 | 0.02 |
| `n2` | 2,015 | 0.02 |
| `root` | 66 | 0.01 |

log-Ne random-walk step scale tau = 1.565

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -768.1 | 222 | 20.68 |
| SNP | +25.2 | 6 | 6.12 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.11 | +0.41 | -0.62 |
| **IBS** | +0.41 | +0.10 | -0.94 |
| **TSI** | -0.62 | -0.94 | +2.07 |

![spectrum](spectrum_fit.png)
