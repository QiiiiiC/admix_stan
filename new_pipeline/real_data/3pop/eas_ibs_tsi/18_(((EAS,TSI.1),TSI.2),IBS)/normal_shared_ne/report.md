# `(((EAS,TSI.1),TSI.2),IBS)`

**Normal, shared Ne** | topology 18 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -39198.71 | +- 0.81 (MC) |
| logZ (importance sampling) | -39152.45 | |
| ESS of the IS weights | 3.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 1 | 18 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 188.5 +- 0.4 | 188.5 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 1.8 +- 0.1 | 190.3 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 1.0 +- 0.0 | 191.3 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 192.3 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `TSI.1`; 1.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,490,624 | 0.06 |
| `IBS` | 229,173 | 0.08 |
| `TSI` | 601,992 | 0.10 |
| `TSI.1` | 5,571 | 0.07 |
| `TSI.2` | 83 | 0.07 |
| `n1` | 2,623 | 0.12 |
| `n2` | 4,961 | 0.06 |
| `root` | 10,063 | 0.02 |

log-Ne random-walk step scale tau = 2.859

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -2,867.6 | 222 | 64.60 |
| SNP | -36,111.5 | 6 | 12051.71 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +113.67 | -122.06 | -102.54 |
| **IBS** | -122.06 | +96.56 | +145.65 |
| **TSI** | -102.54 | +145.65 | +58.49 |

![spectrum](spectrum_fit.png)
