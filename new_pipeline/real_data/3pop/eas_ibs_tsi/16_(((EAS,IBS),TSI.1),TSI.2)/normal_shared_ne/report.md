# `(((EAS,IBS),TSI.1),TSI.2)`

**Normal, shared Ne** | topology 16 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -39192.25 | +- 0.44 (MC) |
| logZ (importance sampling) | -39166.59 | |
| ESS of the IS weights | 3.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 7 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 187.2 +- 0.5 | 187.2 |
| 2 | MERGE | EAS + IBS -> n1 | 1.1 +- 0.0 | 188.3 |
| 3 | MERGE | TSI.2 + n1 -> n2 | 1.1 +- 0.0 | 189.4 |
| 4 | MERGE | TSI.1 + n2 -> root | 1.3 +- 0.0 | 190.6 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,518,248 | 0.03 |
| `IBS` | 237,995 | 0.03 |
| `TSI` | 596,940 | 0.06 |
| `TSI.1` | 161 | 0.09 |
| `TSI.2` | 18,485 | 0.01 |
| `n1` | 1,680 | 0.01 |
| `n2` | 18,605 | 0.01 |
| `root` | 11,130 | 0.01 |

log-Ne random-walk step scale tau = 3.107

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -2,696.8 | 222 | 63.05 |
| SNP | -36,206.3 | 6 | 12083.31 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +113.75 | -122.52 | -102.24 |
| **IBS** | -122.52 | +98.05 | +145.00 |
| **TSI** | -102.24 | +145.00 | +58.54 |

![spectrum](spectrum_fit.png)
