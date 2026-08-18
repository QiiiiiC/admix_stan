# `((CEU,TSI.1),(IBS,TSI.2))`

Model 15 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -601.61 | +- 0.33 (MC) |
| logZ (importance sampling) | -581.94 | |
| ESS of the IS weights | 1.4 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 1 | 22 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 15.2 +- 1.9 | 15.2 |
| 2 | MERGE | TSI.1 + IBS -> n1 | 186.6 +- 8.0 | 201.8 |
| 3 | MERGE | TSI.2 + CEU -> n2 | 874.7 +- 23.7 | 1,076.5 |
| 4 | MERGE | n1 + n2 -> root | 1.3 +- 0.0 | 1,077.8 |

## Admixture fraction

**f = 0.002 +- 0.000** (fraction from `TSI.1`; 0.998 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 605,624 | 0.17 |
| `TSI` | 1,233,664 | 0.14 |
| `CEU` | 303,443 | 0.02 |
| `TSI.1` | 17,950,631 | 0.33 |
| `TSI.2` | 377,093 | 0.05 |
| `n1` | 17,763,715 | 0.33 |
| `n2` | 34,161 | 0.03 |
| `root` | 31,554 | 0.02 |

log-Ne random-walk step scale tau = 1.215

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -528.9 | 110 | 36.33 |
| SNP | +40.5 | 6 | 5.23 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | -0.54 | +2.68 | -1.40 |
| **TSI** | +2.68 | -0.85 | -0.29 |
| **CEU** | -1.40 | -0.29 | +0.93 |

![spectrum](spectrum_fit.png)
