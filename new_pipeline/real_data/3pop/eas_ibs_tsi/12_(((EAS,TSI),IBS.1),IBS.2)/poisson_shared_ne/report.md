# `(((EAS,TSI),IBS.1),IBS.2)`

**Poisson, shared Ne** | topology 12 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -43490.94 | +- 0.23 (MC) |
| logZ (importance sampling) | -43475.74 | |
| ESS of the IS weights | 5.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 7 | 18 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 212.2 +- 1.1 | 212.2 |
| 2 | MERGE | EAS + TSI -> n1 | 1.0 +- 0.0 | 213.2 |
| 3 | MERGE | IBS.2 + n1 -> n2 | 1.0 +- 0.0 | 214.2 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.0 +- 0.0 | 215.2 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `IBS.1`; 1.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,184,206 | 0.02 |
| `IBS` | 438,370 | 0.11 |
| `TSI` | 333,398 | 0.03 |
| `IBS.1` | 16,596 | 0.02 |
| `IBS.2` | 73 | 0.02 |
| `n1` | 642 | 0.02 |
| `n2` | 1,514 | 0.02 |
| `root` | 15,107 | 0.02 |

log-Ne random-walk step scale tau = 3.389

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,124.4 | 222 | 110.40 |
| SNP | -36,105.0 | 6 | 12049.54 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +112.63 | -96.95 | -125.77 |
| **IBS** | -96.95 | +41.27 | +153.67 |
| **TSI** | -125.77 | +153.67 | +95.25 |

![spectrum](spectrum_fit.png)
