# `(((CEU.1,IBS),CEU.2),TSI)`

Model 16 of 21 | admixed leaf: **CEU** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -652.98 | +- 0.27 (MC) |
| logZ (importance sampling) | -632.79 | |
| ESS of the IS weights | 2.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 23 | 21 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | CEU -> CEU.1 + CEU.2 | 143.2 +- 0.8 | 143.2 |
| 2 | MERGE | CEU.2 + IBS -> n1 | 5.5 +- 0.6 | 148.8 |
| 3 | MERGE | CEU.1 + n1 -> n2 | 1.3 +- 0.1 | 150.0 |
| 4 | MERGE | n2 + TSI -> root | 7.0 +- 0.7 | 157.0 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `CEU.1`; 0.000 from `CEU.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 494,336 | 0.07 |
| `TSI` | 420,588 | 0.05 |
| `CEU` | 691,681 | 0.12 |
| `CEU.1` | 2,085 | 0.11 |
| `CEU.2` | 3,797 | 0.09 |
| `n1` | 3,765 | 0.10 |
| `n2` | 3,417 | 0.09 |
| `root` | 13,669 | 0.03 |

log-Ne random-walk step scale tau = 1.629

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -564.2 | 222 | 3.18 |
| SNP | +37.7 | 6 | 6.18 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | -2.10 | +0.96 | +1.21 |
| **TSI** | +0.96 | +2.64 | -3.11 |
| **CEU** | +1.21 | -3.11 | +1.92 |

![spectrum](spectrum_fit.png)
