# `(((EAS,IBS.1),TSI),IBS.2)`

**Poisson, shared Ne** | topology 11 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -4445.75 | +- 0.21 (MC) |
| logZ (importance sampling) | -4428.36 | |
| ESS of the IS weights | 1.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| mode kept / MAP start / runtime | 1 / 11 | 18 s |
| mode search | 10/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 64.7 +- 1.1 | 64.7 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 224.7 +- 0.5 | 289.4 |
| 3 | MERGE | n1 + TSI -> n2 | 13.6 +- 0.0 | 302.9 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.0 +- 0.0 | 304.0 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `IBS.1`; 0.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,725,630 | 0.02 |
| `IBS` | 1,151,820 | 0.05 |
| `TSI` | 322,575 | 0.02 |
| `IBS.1` | 75,905 | 0.04 |
| `IBS.2` | 71 | 0.01 |
| `n1` | 76 | 0.01 |
| `n2` | 1,592 | 0.01 |
| `root` | 784 | 0.01 |

log-Ne random-walk step scale tau = 2.019

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,284.8 | 222 | 23958.62 |
| SNP | +17.8 | 6 | 8.61 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.61 | +0.85 | -2.07 |
| **IBS** | +0.85 | -5.05 | +3.68 |
| **TSI** | -2.07 | +3.68 | +0.49 |

![spectrum](spectrum_fit.png)
