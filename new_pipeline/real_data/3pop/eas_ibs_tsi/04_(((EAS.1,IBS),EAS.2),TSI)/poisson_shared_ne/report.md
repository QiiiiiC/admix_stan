# `(((EAS.1,IBS),EAS.2),TSI)`

**Poisson, shared Ne** | topology 04 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5037.02 | +- 0.24 (MC) |
| logZ (importance sampling) | -5021.36 | |
| ESS of the IS weights | 2.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| mode kept / MAP start / runtime | 3 / 6 | 16 s |
| mode search | 10/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 289.5 +- 0.4 | 289.5 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 12.1 +- 0.1 | 301.6 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 1.2 +- 0.0 | 302.8 |
| 4 | MERGE | n2 + TSI -> root | 1.1 +- 0.0 | 303.9 |

## Admixture fraction

**f = 0.968 +- 0.001** (fraction from `EAS.1`; 0.032 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,724,963 | 0.03 |
| `IBS` | 237,491 | 0.03 |
| `TSI` | 322,219 | 0.03 |
| `EAS.1` | 70 | 0.01 |
| `EAS.2` | 1,718 | 0.01 |
| `n1` | 1,679 | 0.01 |
| `n2` | 1,299 | 0.01 |
| `root` | 768 | 0.01 |

log-Ne random-walk step scale tau = 1.677

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,895.6 | 222 | 23235.13 |
| SNP | +33.2 | 6 | 3.46 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.45 | -0.29 | -0.61 |
| **IBS** | -0.29 | -1.30 | +1.99 |
| **TSI** | -0.61 | +1.99 | -0.72 |

![spectrum](spectrum_fit.png)
