# `((CEU.1,IBS),(CEU.2,TSI))`

Model 21 of 21 | admixed leaf: **CEU** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +1284.07 | +- 0.15 (MC) |
| logZ (importance sampling) | +1300.12 | |
| ESS of the IS weights | 1.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 7 | 22 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | CEU -> CEU.1 + CEU.2 | 127.6 +- 0.7 | 127.6 |
| 2 | MERGE | CEU.1 + IBS -> n1 | 2.2 +- 0.1 | 129.8 |
| 3 | MERGE | CEU.2 + TSI -> n2 | 31.3 +- 0.7 | 161.1 |
| 4 | MERGE | n1 + n2 -> root | 5.0 +- 0.2 | 166.1 |

## Admixture fraction

**f = 0.880 +- 0.001** (fraction from `CEU.1`; 0.120 from `CEU.2`)

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 744,744 | 0.08 |
| `TSI` | 506,019 | 0.08 |
| `CEU` | 843,140 | 0.19 |
| `CEU.1` | 22,922 | 0.07 |
| `CEU.2` | 131 | 0.06 |
| `n1` | 19,826 | 0.04 |
| `n2` | 6,361 | 0.05 |
| `root` | 10,613 | 0.03 |

log-Ne random-walk step scale tau = 1.700

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +1,350.7 | 110 | 2.15 |
| SNP | +50.7 | 6 | 1.85 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | -0.99 | +0.67 | +0.42 |
| **TSI** | +0.67 | +0.72 | -1.02 |
| **CEU** | +0.42 | -1.02 | +0.62 |

![spectrum](spectrum_fit.png)
