# `(((IBS.1,TSI),CEU),IBS.2)`

Model 07 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -2140.58 | +- 0.13 (MC) |
| logZ (importance sampling) | -2130.14 | |
| ESS of the IS weights | 1.8 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 47 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 42.2 +- 0.9 | 42.2 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 161.4 +- 0.9 | 203.7 |
| 3 | MERGE | n1 + CEU -> n2 | 1.2 +- 0.0 | 204.9 |
| 4 | MERGE | IBS.1 + n2 -> root | 2.1 +- 0.0 | 207.0 |

## Admixture fraction

**f = 0.443 +- 0.008** (fraction from `IBS.1`; 0.557 from `IBS.2`)

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 1,351,644 | 0.07 |
| `TSI` | 497,382 | 0.04 |
| `CEU` | 276,633 | 0.03 |
| `IBS.1` | 102,451 | 0.05 |
| `IBS.2` | 85,686 | 0.05 |
| `n1` | 285 | 0.04 |
| `n2` | 501 | 0.04 |
| `root` | 2,202 | 0.02 |

log-Ne random-walk step scale tau = 2.325

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -2,007.6 | 222 | 25.75 |
| SNP | +21.2 | 6 | 11.68 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | -0.84 | +1.98 | -0.63 |
| **TSI** | +1.98 | +4.10 | -5.04 |
| **CEU** | -0.63 | -5.04 | +4.39 |

![spectrum](spectrum_fit.png)
