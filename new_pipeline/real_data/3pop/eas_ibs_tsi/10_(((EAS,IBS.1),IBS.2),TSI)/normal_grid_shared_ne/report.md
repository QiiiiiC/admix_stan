# `(((EAS,IBS.1),IBS.2),TSI)`

**Normal, shared Ne, recent grid** | topology 10 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -35373.21 | +- 1.38 (MC) |
| logZ (importance sampling) | -35295.26 | |
| ESS of the IS weights | 2.4 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 1 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 92.0 +- 4.4 | 92.0 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 75.3 +- 6.0 | 167.3 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 1.5 +- 0.1 | 168.8 |
| 4 | MERGE | n2 + TSI -> root | 1.0 +- 0.0 | 169.8 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `IBS.1`; 0.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,627,729 | 9,530,988 |
| `IBS` | 16,938,313 | 7,920,381 |
| `TSI` | 3,299,554 | 2,525,955 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,402,747 | 0.08 |
| `IBS` | 664,007 | 0.10 |
| `TSI` | 313,610 | 0.05 |
| `IBS.1` | 39,049 | 0.16 |
| `IBS.2` | 0 | 0.24 |
| `n1` | 827 | 0.03 |
| `n2` | 52,586 | 0.03 |
| `root` | 35,653 | 0.07 |

log-Ne random-walk step scale tau = 2.615

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +396.8 | 222 | 35.34 |
| SNP | -35,436.5 | 6 | 11826.72 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +111.55 | -95.97 | -124.62 |
| **IBS** | -95.97 | +40.93 | +152.02 |
| **TSI** | -124.62 | +152.02 | +94.59 |

![spectrum](spectrum_fit.png)
