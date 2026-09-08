# `(((IBS,TSI.1),EAS),TSI.2)`

**Poisson, shared Ne** | topology 19 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5881.29 | +- 0.17 (MC) |
| logZ (importance sampling) | -5867.03 | |
| ESS of the IS weights | 1.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| mode kept / MAP start / runtime | 3 / 9 | 19 s |
| mode search | 10/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 167.4 +- 0.7 | 167.4 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 69.5 +- 0.2 | 237.0 |
| 3 | MERGE | n1 + EAS -> n2 | 160.9 +- 0.3 | 397.9 |
| 4 | MERGE | TSI.1 + n2 -> root | 344.3 +- 0.5 | 742.1 |

## Admixture fraction

**f = 0.007 +- 0.000** (fraction from `TSI.1`; 0.993 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 898,053 | 0.02 |
| `IBS` | 315,230 | 0.02 |
| `TSI` | 420,998 | 0.05 |
| `TSI.1` | 27,412 | 0.00 |
| `TSI.2` | 254,653 | 0.03 |
| `n1` | 891 | 0.01 |
| `n2` | 8 | 0.01 |
| `root` | 9,285 | 0.00 |

log-Ne random-walk step scale tau = 1.806

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,760.3 | 222 | 91.80 |
| SNP | +38.8 | 6 | 1.60 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.21 | +0.18 | -0.58 |
| **IBS** | +0.18 | -1.27 | +1.00 |
| **TSI** | -0.58 | +1.00 | +0.17 |

![spectrum](spectrum_fit.png)
