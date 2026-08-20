# `(((CEU,TSI),IBS.1),IBS.2)`

Model 06 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -2698.44 | +- 0.23 (MC) |
| logZ (importance sampling) | -2682.99 | |
| ESS of the IS weights | 7.9 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 13 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 107.6 +- 0.9 | 107.6 |
| 2 | MERGE | TSI + CEU -> n1 | 44.7 +- 1.0 | 152.3 |
| 3 | MERGE | IBS.2 + n1 -> n2 | 1.3 +- 0.0 | 153.6 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.6 +- 0.0 | 155.2 |

## Admixture fraction

**f = 0.156 +- 0.028** (fraction from `IBS.1`; 0.844 from `IBS.2`)

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 647,281 | 0.14 |
| `TSI` | 434,880 | 0.06 |
| `CEU` | 245,084 | 0.04 |
| `IBS.1` | 10,101 | 0.06 |
| `IBS.2` | 24,707 | 0.05 |
| `n1` | 13,020 | 0.04 |
| `n2` | 12,081 | 0.04 |
| `root` | 11,378 | 0.04 |

log-Ne random-walk step scale tau = 1.109

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,451.2 | 222 | 10.76 |
| SNP | -1,143.0 | 6 | 399.72 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | +0.58 | +5.85 | -4.62 |
| **TSI** | +5.85 | +26.94 | -30.02 |
| **CEU** | -4.62 | -30.02 | +26.62 |

![spectrum](spectrum_fit.png)
