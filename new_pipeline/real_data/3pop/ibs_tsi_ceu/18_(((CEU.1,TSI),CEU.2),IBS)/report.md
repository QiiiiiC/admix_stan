# `(((CEU.1,TSI),CEU.2),IBS)`

Model 18 of 21 | admixed leaf: **CEU** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1467.58 | +- 0.27 (MC) |
| logZ (importance sampling) | -1450.63 | |
| ESS of the IS weights | 3.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 13 | 19 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | CEU -> CEU.1 + CEU.2 | 122.0 +- 1.0 | 122.0 |
| 2 | MERGE | CEU.2 + TSI -> n1 | 38.7 +- 1.7 | 160.7 |
| 3 | MERGE | CEU.1 + n1 -> n2 | 5.3 +- 0.2 | 166.0 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 167.1 |

## Admixture fraction

**f = 0.996 +- 0.001** (fraction from `CEU.1`; 0.004 from `CEU.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 349,382 | 0.03 |
| `TSI` | 536,248 | 0.08 |
| `CEU` | 663,296 | 0.07 |
| `CEU.1` | 12,824 | 0.05 |
| `CEU.2` | 2,946 | 0.08 |
| `n1` | 2,835 | 0.05 |
| `n2` | 8,401 | 0.06 |
| `root` | 7,337 | 0.04 |

log-Ne random-walk step scale tau = 1.558

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,362.9 | 222 | 11.32 |
| SNP | +23.2 | 6 | 10.99 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | -0.96 | -2.94 | +2.93 |
| **TSI** | -2.94 | +4.92 | -3.75 |
| **CEU** | +2.93 | -3.75 | +1.58 |

![spectrum](spectrum_fit.png)
