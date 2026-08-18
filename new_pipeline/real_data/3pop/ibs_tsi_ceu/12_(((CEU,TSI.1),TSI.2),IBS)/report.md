# `(((CEU,TSI.1),TSI.2),IBS)`

Model 12 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +1025.78 | +- 0.13 (MC) |
| logZ (importance sampling) | +1034.90 | |
| ESS of the IS weights | 9.3 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 23 | 21 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 143.3 +- 0.4 | 143.3 |
| 2 | MERGE | TSI.2 + CEU -> n1 | 1.5 +- 0.1 | 144.8 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 14.6 +- 0.3 | 159.4 |
| 4 | MERGE | n2 + IBS -> root | 1.6 +- 0.2 | 161.0 |

## Admixture fraction

**f = 0.980 +- 0.004** (fraction from `TSI.1`; 0.020 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 643,072 | 0.10 |
| `TSI` | 791,640 | 0.08 |
| `CEU` | 573,599 | 0.08 |
| `TSI.1` | 6,627 | 0.06 |
| `TSI.2` | 4,082 | 0.03 |
| `n1` | 4,084 | 0.03 |
| `n2` | 9,915 | 0.03 |
| `root` | 10,828 | 0.03 |

log-Ne random-walk step scale tau = 1.474

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +1,087.0 | 110 | 7.00 |
| SNP | +43.5 | 6 | 4.25 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | -0.80 | -2.11 | +2.20 |
| **TSI** | -2.11 | +2.02 | -1.16 |
| **CEU** | +2.20 | -1.16 | -0.16 |

![spectrum](spectrum_fit.png)
