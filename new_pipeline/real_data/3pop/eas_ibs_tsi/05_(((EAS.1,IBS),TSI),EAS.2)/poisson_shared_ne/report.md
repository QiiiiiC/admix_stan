# `(((EAS.1,IBS),TSI),EAS.2)`

**Poisson, shared Ne** | topology 05 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -6163.62 | +- 3.64 (MC) |
| logZ (importance sampling) | -6001.41 | |
| ESS of the IS weights | 1.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 1 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 260.4 +- 1.9 | 261.4 |
| 3 | MERGE | n1 + TSI -> n2 | 1.0 +- 0.0 | 262.4 |
| 4 | MERGE | EAS.1 + n2 -> root | 51.7 +- 2.0 | 314.1 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 695,084 | 0.07 |
| `IBS` | 325,912 | 0.21 |
| `TSI` | 429,336 | 0.22 |
| `EAS.1` | 888,380 | 0.06 |
| `EAS.2` | 34 | 0.23 |
| `n1` | 319 | 0.04 |
| `n2` | 290 | 0.04 |
| `root` | 894 | 0.06 |

log-Ne random-walk step scale tau = 1.649

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,953.3 | 222 | 219.35 |
| SNP | +6.9 | 6 | 12.24 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.61 | +1.06 | -2.28 |
| **IBS** | +1.06 | -5.20 | +3.42 |
| **TSI** | -2.28 | +3.42 | +1.14 |

![spectrum](spectrum_fit.png)
