# `(((CEU,IBS.1),TSI),IBS.2)`

Model 05 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -845.73 | +- 2.09 (MC) |
| logZ (importance sampling) | -739.55 | |
| ESS of the IS weights | 2.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 23 | 24 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | IBS.2 + CEU -> n1 | 1.2 +- 0.0 | 2.2 |
| 3 | MERGE | n1 + TSI -> n2 | 164.4 +- 4.0 | 166.6 |
| 4 | MERGE | IBS.1 + n2 -> root | 137.3 +- 4.6 | 303.9 |

## Admixture fraction

**f = 0.944 +- 0.003** (fraction from `IBS.1`; 0.056 from `IBS.2`)

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 385,889 | 0.23 |
| `TSI` | 303,341 | 0.13 |
| `CEU` | 274,059 | 0.13 |
| `IBS.1` | 409,420 | 0.23 |
| `IBS.2` | 282,564 | 0.13 |
| `n1` | 279,824 | 0.13 |
| `n2` | 288,208 | 0.13 |
| `root` | 50 | 0.16 |

log-Ne random-walk step scale tau = 1.239

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +550.2 | 110 | 15.47 |
| SNP | -1,191.8 | 6 | 416.00 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | +10.83 | +1.78 | -10.97 |
| **TSI** | +1.78 | +26.24 | -27.56 |
| **CEU** | -10.97 | -27.56 | +27.79 |

![spectrum](spectrum_fit.png)
