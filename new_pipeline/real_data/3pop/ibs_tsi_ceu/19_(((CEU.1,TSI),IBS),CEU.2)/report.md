# `(((CEU.1,TSI),IBS),CEU.2)`

Model 19 of 21 | admixed leaf: **CEU** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1478.82 | +- 0.08 (MC) |
| logZ (importance sampling) | -1471.40 | |
| ESS of the IS weights | 7.4 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 13 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | CEU -> CEU.1 + CEU.2 | 120.0 +- 0.5 | 120.0 |
| 2 | MERGE | CEU.2 + TSI -> n1 | 41.5 +- 0.2 | 161.5 |
| 3 | MERGE | n1 + IBS -> n2 | 5.0 +- 0.1 | 166.5 |
| 4 | MERGE | CEU.1 + n2 -> root | 1.0 +- 0.0 | 167.5 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `CEU.1`; 0.000 from `CEU.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 349,244 | 0.04 |
| `TSI` | 537,080 | 0.04 |
| `CEU` | 661,110 | 0.04 |
| `CEU.1` | 14,058 | 0.01 |
| `CEU.2` | 2,556 | 0.06 |
| `n1` | 2,531 | 0.02 |
| `n2` | 11,805 | 0.02 |
| `root` | 7,206 | 0.01 |

log-Ne random-walk step scale tau = 1.869

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,375.6 | 222 | 11.40 |
| SNP | +34.4 | 6 | 7.28 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | +0.65 | -3.62 | +1.96 |
| **TSI** | -3.62 | +4.05 | -2.57 |
| **CEU** | +1.96 | -2.57 | +1.11 |

![spectrum](spectrum_fit.png)
