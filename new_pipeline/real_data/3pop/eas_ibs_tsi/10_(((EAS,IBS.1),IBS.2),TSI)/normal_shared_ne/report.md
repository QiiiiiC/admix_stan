# `(((EAS,IBS.1),IBS.2),TSI)`

**Normal, shared Ne** | topology 10 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -38276.81 | +- 0.14 (MC) |
| logZ (importance sampling) | -38265.28 | |
| ESS of the IS weights | 1.4 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 1 | 18 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 183.9 +- 1.1 | 183.9 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 2.3 +- 0.0 | 186.2 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 1.0 +- 0.0 | 187.2 |
| 4 | MERGE | n2 + TSI -> root | 1.1 +- 0.0 | 188.3 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `IBS.1`; 1.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,517,759 | 0.03 |
| `IBS` | 719,555 | 0.04 |
| `TSI` | 340,557 | 0.02 |
| `IBS.1` | 7,204 | 0.03 |
| `IBS.2` | 100 | 0.02 |
| `n1` | 2,103 | 0.02 |
| `n2` | 4,935 | 0.03 |
| `root` | 12,268 | 0.03 |

log-Ne random-walk step scale tau = 2.806

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,904.5 | 222 | 55.97 |
| SNP | -36,177.3 | 6 | 12073.65 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +113.54 | -101.01 | -123.48 |
| **IBS** | -101.01 | +55.50 | +146.70 |
| **TSI** | -123.48 | +146.70 | +97.43 |

![spectrum](spectrum_fit.png)
