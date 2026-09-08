# `(((EAS,IBS.1),TSI),IBS.2)`

**Normal, shared Ne** | topology 11 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -35440.67 | +- 0.41 (MC) |
| logZ (importance sampling) | -35409.23 | |
| ESS of the IS weights | 1.3 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| mode kept / MAP start / runtime | 1 / 2 | 20 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 128.1 +- 1.0 | 128.1 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 28.9 +- 1.2 | 157.1 |
| 3 | MERGE | n1 + TSI -> n2 | 6.6 +- 0.1 | 163.7 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.0 +- 0.0 | 164.7 |

## Admixture fraction

**f = 0.978 +- 0.002** (fraction from `IBS.1`; 0.022 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,577,393 | 0.04 |
| `IBS` | 596,402 | 0.06 |
| `TSI` | 328,108 | 0.02 |
| `IBS.1` | 12,729 | 0.05 |
| `IBS.2` | 1 | 0.16 |
| `n1` | 4,506 | 0.03 |
| `n2` | 142,825 | 0.01 |
| `root` | 43,304 | 0.01 |

log-Ne random-walk step scale tau = 1.442

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +436.8 | 222 | 34.71 |
| SNP | -35,598.1 | 6 | 11880.58 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +111.86 | -96.34 | -124.85 |
| **IBS** | -96.34 | +41.57 | +152.11 |
| **TSI** | -124.85 | +152.11 | +94.96 |

![spectrum](spectrum_fit.png)
