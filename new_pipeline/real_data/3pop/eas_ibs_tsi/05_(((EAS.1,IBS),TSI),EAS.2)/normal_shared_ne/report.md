# `(((EAS.1,IBS),TSI),EAS.2)`

**Normal, shared Ne** | topology 05 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +4023.16 | +- 0.13 (MC) |
| logZ (importance sampling) | +4031.04 | |
| ESS of the IS weights | 20.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| mode kept / MAP start / runtime | 1 / 6 | 19 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 127.4 +- 0.6 | 127.4 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 3.8 +- 0.1 | 131.2 |
| 3 | MERGE | n1 + TSI -> n2 | 19.6 +- 1.2 | 150.8 |
| 4 | MERGE | EAS.1 + n2 -> root | 2,546.0 +- 30.9 | 2,696.8 |

## Admixture fraction

**f = 0.993 +- 0.000** (fraction from `EAS.1`; 0.007 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,583,574 | 0.03 |
| `IBS` | 550,671 | 0.03 |
| `TSI` | 408,925 | 0.03 |
| `EAS.1` | 44,096 | 0.02 |
| `EAS.2` | 11,848 | 0.10 |
| `n1` | 9,166 | 0.06 |
| `n2` | 20,419 | 0.02 |
| `root` | 14,517 | 0.09 |

log-Ne random-walk step scale tau = 1.339

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +4,128.5 | 222 | 1.48 |
| SNP | +36.3 | 6 | 2.44 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.04 | +0.97 | -1.06 |
| **IBS** | +0.97 | -1.93 | +0.09 |
| **TSI** | -1.06 | +0.09 | +1.93 |

![spectrum](spectrum_fit.png)
