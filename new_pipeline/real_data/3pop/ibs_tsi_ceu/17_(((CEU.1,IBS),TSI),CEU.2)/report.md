# `(((CEU.1,IBS),TSI),CEU.2)`

Model 17 of 21 | admixed leaf: **CEU** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -364.39 | +- 1.23 (MC) |
| logZ (importance sampling) | -304.60 | |
| ESS of the IS weights | 2.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 7 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | CEU -> CEU.1 + CEU.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | CEU.2 + IBS -> n1 | 2.2 +- 0.6 | 3.2 |
| 3 | MERGE | n1 + TSI -> n2 | 297.6 +- 0.9 | 300.8 |
| 4 | MERGE | CEU.1 + n2 -> root | 1.1 +- 0.0 | 301.9 |

## Admixture fraction

**f = 0.011 +- 0.000** (fraction from `CEU.1`; 0.989 from `CEU.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 4,734,563 | 0.25 |
| `TSI` | 263,570 | 0.09 |
| `CEU` | 4,128,929 | 0.25 |
| `CEU.1` | 10 | 0.17 |
| `CEU.2` | 4,734,072 | 0.24 |
| `n1` | 4,680,410 | 0.24 |
| `n2` | 50 | 0.21 |
| `root` | 46 | 0.04 |

log-Ne random-walk step scale tau = 1.611

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +288.1 | 110 | 20.25 |
| SNP | -462.3 | 6 | 172.83 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | +13.70 | -21.47 | +2.77 |
| **TSI** | -21.47 | +15.34 | -6.42 |
| **CEU** | +2.77 | -6.42 | +3.82 |

![spectrum](spectrum_fit.png)
