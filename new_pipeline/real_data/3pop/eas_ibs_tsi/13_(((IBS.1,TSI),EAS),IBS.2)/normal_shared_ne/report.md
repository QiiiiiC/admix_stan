# `(((IBS.1,TSI),EAS),IBS.2)`

**Normal, shared Ne** | topology 13 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1584.73 | +- 0.76 (MC) |
| logZ (importance sampling) | -1539.31 | |
| ESS of the IS weights | 11.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| mode kept / MAP start / runtime | 1 / 1 | 20 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 31.5 +- 0.7 | 31.5 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 181.0 +- 0.8 | 212.5 |
| 3 | MERGE | n1 + EAS -> n2 | 198.2 +- 4.3 | 410.7 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.2 +- 0.0 | 411.9 |

## Admixture fraction

**f = 0.003 +- 0.000** (fraction from `IBS.1`; 0.997 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 757,884 | 0.03 |
| `IBS` | 3,371,885 | 0.04 |
| `TSI` | 483,178 | 0.05 |
| `IBS.1` | 1 | 0.05 |
| `IBS.2` | 245,797 | 0.03 |
| `n1` | 1,119 | 0.03 |
| `n2` | 2 | 0.03 |
| `root` | 2 | 0.03 |

log-Ne random-walk step scale tau = 1.635

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,387.6 | 222 | 52.34 |
| SNP | +8.5 | 6 | 11.72 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +1.97 | -1.20 | -2.70 |
| **IBS** | -1.20 | -0.24 | +2.71 |
| **TSI** | -2.70 | +2.71 | +2.60 |

![spectrum](spectrum_fit.png)
