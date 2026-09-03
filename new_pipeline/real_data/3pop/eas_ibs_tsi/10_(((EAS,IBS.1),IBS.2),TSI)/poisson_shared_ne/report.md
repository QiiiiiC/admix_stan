# `(((EAS,IBS.1),IBS.2),TSI)`

**Poisson, shared Ne** | topology 10 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5132.59 | +- 0.50 (MC) |
| logZ (importance sampling) | -5104.22 | |
| ESS of the IS weights | 4.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 1 | 15 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 288.7 +- 0.3 | 288.7 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 1.0 +- 0.0 | 289.7 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 13.1 +- 0.0 | 302.8 |
| 4 | MERGE | n2 + TSI -> root | 1.0 +- 0.0 | 303.8 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `IBS.1`; 0.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,726,629 | 0.05 |
| `IBS` | 237,791 | 0.05 |
| `TSI` | 322,382 | 0.03 |
| `IBS.1` | 23,631 | 0.03 |
| `IBS.2` | 73 | 0.02 |
| `n1` | 73 | 0.02 |
| `n2` | 1,029 | 0.02 |
| `root` | 772 | 0.02 |

log-Ne random-walk step scale tau = 1.663

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,908.1 | 222 | 23230.65 |
| SNP | +30.1 | 6 | 4.49 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.43 | -0.41 | -0.44 |
| **IBS** | -0.41 | -1.01 | +1.92 |
| **TSI** | -0.44 | +1.92 | -0.97 |

![spectrum](spectrum_fit.png)
