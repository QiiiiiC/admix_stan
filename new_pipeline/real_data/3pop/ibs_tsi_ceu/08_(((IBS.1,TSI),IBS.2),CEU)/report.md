# `(((IBS.1,TSI),IBS.2),CEU)`

Model 08 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -2151.75 | +- 0.16 (MC) |
| logZ (importance sampling) | -2140.75 | |
| ESS of the IS weights | 2.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 7 | 21 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 40.1 +- 0.8 | 40.1 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 160.5 +- 0.4 | 200.6 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 3.0 +- 0.1 | 203.6 |
| 4 | MERGE | n2 + CEU -> root | 2.4 +- 0.1 | 206.0 |

## Admixture fraction

**f = 0.881 +- 0.005** (fraction from `IBS.1`; 0.119 from `IBS.2`)

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 1,348,468 | 0.09 |
| `TSI` | 492,286 | 0.03 |
| `CEU` | 265,537 | 0.04 |
| `IBS.1` | 470,782 | 0.05 |
| `IBS.2` | 4,207 | 0.09 |
| `n1` | 1,276 | 0.03 |
| `n2` | 1,164 | 0.02 |
| `root` | 1,837 | 0.02 |

log-Ne random-walk step scale tau = 1.557

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -2,039.5 | 222 | 27.18 |
| SNP | +15.0 | 6 | 13.75 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | +0.16 | +3.14 | -2.34 |
| **TSI** | +3.14 | +3.37 | -4.79 |
| **CEU** | -2.34 | -4.79 | +5.04 |

![spectrum](spectrum_fit.png)
