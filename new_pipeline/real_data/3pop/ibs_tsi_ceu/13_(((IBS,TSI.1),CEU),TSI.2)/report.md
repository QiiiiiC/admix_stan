# `(((IBS,TSI.1),CEU),TSI.2)`

Model 13 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -346.96 | +- 0.42 (MC) |
| logZ (importance sampling) | -324.56 | |
| ESS of the IS weights | 8.9 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 7 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 1.1 +- 0.0 | 2.1 |
| 3 | MERGE | n1 + CEU -> n2 | 324.0 +- 0.5 | 326.0 |
| 4 | MERGE | TSI.1 + n2 -> root | 1.6 +- 0.0 | 327.6 |

## Admixture fraction

**f = 0.061 +- 0.003** (fraction from `TSI.1`; 0.939 from `TSI.2`)

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 45,933,606 | 0.34 |
| `TSI` | 20,534,787 | 0.31 |
| `CEU` | 234,976 | 0.07 |
| `TSI.1` | 1,890 | 0.32 |
| `TSI.2` | 33,539,425 | 0.34 |
| `n1` | 31,523,818 | 0.34 |
| `n2` | 3 | 0.04 |
| `root` | 1 | 0.05 |

log-Ne random-walk step scale tau = 1.800

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +490.8 | 110 | 17.64 |
| SNP | -624.5 | 6 | 226.91 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | +16.42 | +10.83 | -22.34 |
| **TSI** | +10.83 | +4.64 | -9.40 |
| **CEU** | -22.34 | -9.40 | +18.73 |

![spectrum](spectrum_fit.png)
