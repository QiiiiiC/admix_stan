# `(((EAS.1,TSI),IBS),EAS.2)`

**Normal, shared Ne** | topology 07 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5194.06 | +- 0.04 (MC) |
| logZ (importance sampling) | -5188.38 | |
| ESS of the IS weights | 3.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 13 | 18 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 1.0 +- 0.0 | 2.0 |
| 3 | MERGE | n1 + IBS -> n2 | 1.0 +- 0.0 | 3.0 |
| 4 | MERGE | EAS.1 + n2 -> root | 40,150.8 +- 322.3 | 40,153.8 |

## Admixture fraction

**f = 0.994 +- 0.000** (fraction from `EAS.1`; 0.006 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 694,352 | 0.01 |
| `IBS` | 320,102 | 0.01 |
| `TSI` | 320,110 | 0.01 |
| `EAS.1` | 697,299 | 0.01 |
| `EAS.2` | 320,117 | 0.01 |
| `n1` | 320,109 | 0.01 |
| `n2` | 320,099 | 0.01 |
| `root` | 359,026 | 0.01 |

log-Ne random-walk step scale tau = 0.026

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,617.5 | 222 | 72.27 |
| SNP | +20.8 | 6 | 7.59 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.04 | +0.03 | -0.10 |
| **IBS** | +0.03 | +3.61 | -3.92 |
| **TSI** | -0.10 | -3.92 | +3.88 |

![spectrum](spectrum_fit.png)
