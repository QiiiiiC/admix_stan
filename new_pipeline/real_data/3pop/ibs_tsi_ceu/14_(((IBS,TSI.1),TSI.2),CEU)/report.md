# `(((IBS,TSI.1),TSI.2),CEU)`

Model 14 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +190.35 | +- 0.43 (MC) |
| logZ (importance sampling) | +214.40 | |
| ESS of the IS weights | 3.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -25.77 | already applied |
| seed kept / runtime | 31 | 24 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 2.5 +- 0.6 | 3.5 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 304.8 +- 1.8 | 308.3 |
| 4 | MERGE | n2 + CEU -> root | 8.5 +- 0.1 | 316.8 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `IBS` | 435,719 | 0.04 |
| `TSI` | 244,137 | 0.11 |
| `CEU` | 376,940 | 0.07 |
| `TSI.1` | 286,201 | 0.11 |
| `TSI.2` | 447,082 | 0.04 |
| `n1` | 470,454 | 0.05 |
| `n2` | 3,226 | 0.09 |
| `root` | 14 | 0.02 |

log-Ne random-walk step scale tau = 2.383

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +663.5 | 110 | 14.26 |
| SNP | -255.7 | 6 | 103.96 |

chi2/n is the mean squared standardized residual: 1.0 = fits inside its own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- different term counts and different units.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | IBS | TSI | CEU |
|---|---|---|---|
| **IBS** | +3.41 | -14.29 | +6.96 |
| **TSI** | -14.29 | +15.43 | -9.60 |
| **CEU** | +6.96 | -9.60 | +4.31 |

![spectrum](spectrum_fit.png)
