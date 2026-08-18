# `((CEU,IBS.1),(IBS.2,TSI))`

Model 09 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +640.33 | +- 0.39 (MC) |
| logZ (importance sampling) | +666.43 | |
| ESS of the IS weights | 1.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 7 | 22 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | IBS.1 + TSI -> n1 | 149.8 +- 1.8 | 150.8 |
| 3 | MERGE | IBS.2 + CEU -> n2 | 3.8 +- 0.9 | 154.6 |
| 4 | MERGE | n1 + n2 -> root | 10.2 +- 1.2 | 164.8 |

## Admixture fraction

**f = 0.583 +- 0.015** (fraction from `IBS.1`; 0.417 from `IBS.2`)

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 420,904 | 0.12 |
| `TSI` | 397,429 | 0.08 |
| `CEU` | 542,057 | 0.11 |
| `IBS.1` | 1,555,076 | 0.19 |
| `IBS.2` | 75,454 | 0.09 |
| `n1` | 163,843 | 0.07 |
| `n2` | 1,749 | 0.13 |
| `root` | 10,898 | 0.04 |

log-Ne random-walk step scale tau = 1.767

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +1,230.0 | 110 | 4.10 |
| SNP | -458.6 | 6 | 171.62 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | +24.91 | -15.01 | -11.81 |
| **TSI** | -15.01 | +1.98 | +4.44 |
| **CEU** | -11.81 | +4.44 | +2.28 |

![spectrum](spectrum_fit.png)
