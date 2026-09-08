# `(((EAS,IBS.1),TSI),IBS.2)`

**Normal, shared Ne, recent grid** | topology 11 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -35320.76 | +- 1.74 (MC) |
| logZ (importance sampling) | -35244.59 | |
| ESS of the IS weights | 4.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 1 / 2 | 113 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 87.9 +- 2.3 | 87.9 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 78.4 +- 1.9 | 166.3 |
| 3 | MERGE | n1 + TSI -> n2 | 2.4 +- 0.1 | 168.7 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.0 +- 0.0 | 169.7 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `IBS.1`; 0.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 3,553,366 | 4,101,159 |
| `IBS` | 2,380,343 | 2,079,281 |
| `TSI` | 3,426,188 | 2,059,866 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,402,842 | 0.07 |
| `IBS` | 725,046 | 0.19 |
| `TSI` | 314,342 | 0.08 |
| `IBS.1` | 42,838 | 0.07 |
| `IBS.2` | 0 | 0.30 |
| `n1` | 1,298 | 0.06 |
| `n2` | 66,372 | 0.10 |
| `root` | 35,668 | 0.02 |

log-Ne random-walk step scale tau = 2.887

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +413.6 | 222 | 35.17 |
| SNP | -35,460.3 | 6 | 11834.63 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +111.60 | -96.04 | -124.63 |
| **IBS** | -96.04 | +41.18 | +151.91 |
| **TSI** | -124.63 | +151.91 | +94.72 |

![spectrum](spectrum_fit.png)
