# `(((IBS,TSI),CEU.1),CEU.2)`

Model 20 of 21 | admixed leaf: **CEU** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -782.15 | +- 2.01 (MC) |
| logZ (importance sampling) | -716.22 | |
| ESS of the IS weights | 1.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 31 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | CEU -> CEU.1 + CEU.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | IBS + TSI -> n1 | 1.0 +- 0.0 | 2.0 |
| 3 | MERGE | CEU.2 + n1 -> n2 | 144.1 +- 1.9 | 146.1 |
| 4 | MERGE | CEU.1 + n2 -> root | 1.3 +- 0.0 | 147.4 |

## Admixture fraction

**f = 0.003 +- 0.000** (fraction from `CEU.1`; 0.997 from `CEU.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 371,183,660 | 1.17 |
| `TSI` | 370,331,648 | 1.17 |
| `CEU` | 390,907 | 0.17 |
| `CEU.1` | 0 | 0.28 |
| `CEU.2` | 404,000 | 0.17 |
| `n1` | 367,885,611 | 1.17 |
| `n2` | 14,740 | 0.13 |
| `root` | 15,366 | 0.07 |

log-Ne random-walk step scale tau = 1.874

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +650.7 | 110 | 15.19 |
| SNP | -1,263.9 | 6 | 440.05 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | +18.96 | -35.61 | +7.97 |
| **TSI** | -35.61 | +26.96 | -12.18 |
| **CEU** | +7.97 | -12.18 | +5.90 |

![spectrum](spectrum_fit.png)
