# `(((CEU,IBS),TSI.1),TSI.2)`

Model 10 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1928.89 | +- 0.31 (MC) |
| logZ (importance sampling) | -1910.72 | |
| ESS of the IS weights | 6.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 7 | 21 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 13.0 +- 0.9 | 13.0 |
| 2 | MERGE | IBS + CEU -> n1 | 149.2 +- 0.9 | 162.2 |
| 3 | MERGE | TSI.2 + n1 -> n2 | 1.2 +- 0.0 | 163.3 |
| 4 | MERGE | TSI.1 + n2 -> root | 14.8 +- 0.2 | 178.2 |

## Admixture fraction

**f = 0.542 +- 0.007** (fraction from `TSI.1`; 0.458 from `TSI.2`)

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 425,948 | 0.05 |
| `TSI` | 1,498,308 | 0.23 |
| `CEU` | 290,558 | 0.06 |
| `TSI.1` | 6,821,653 | 0.39 |
| `TSI.2` | 73,775 | 0.05 |
| `n1` | 2,901 | 0.06 |
| `n2` | 2,486 | 0.05 |
| `root` | 29,090 | 0.02 |

log-Ne random-walk step scale tau = 1.743

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -819.1 | 222 | 4.93 |
| SNP | -987.5 | 6 | 347.91 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | +14.71 | +16.63 | -24.87 |
| **TSI** | +16.63 | +8.61 | -15.95 |
| **CEU** | -24.87 | -15.95 | +25.30 |

![spectrum](spectrum_fit.png)
