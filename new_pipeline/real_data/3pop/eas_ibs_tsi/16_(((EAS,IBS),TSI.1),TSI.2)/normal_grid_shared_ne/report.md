# `(((EAS,IBS),TSI.1),TSI.2)`

**Normal, shared Ne, recent grid** | topology 16 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -37639.91 | +- 1.12 (MC) |
| logZ (importance sampling) | -37580.74 | |
| ESS of the IS weights | 4.9 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 7 | 24 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 167.3 +- 1.0 | 167.3 |
| 2 | MERGE | EAS + IBS -> n1 | 1.3 +- 0.0 | 168.5 |
| 3 | MERGE | TSI.2 + n1 -> n2 | 1.1 +- 0.0 | 169.7 |
| 4 | MERGE | TSI.1 + n2 -> root | 1.1 +- 0.0 | 170.8 |

## Admixture fraction

**f = 0.041 +- 0.002** (fraction from `TSI.1`; 0.959 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 776,327,228 | 516,145,141 |
| `IBS` | 5,606,839 | 4,010,888 |
| `TSI` | 24,149,479 | 1,910,378 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,243,016 | 0.06 |
| `IBS` | 222,418 | 0.07 |
| `TSI` | 402,899 | 0.03 |
| `TSI.1` | 0 | 0.07 |
| `TSI.2` | 25,982 | 0.02 |
| `n1` | 43,335 | 0.01 |
| `n2` | 3,348 | 0.02 |
| `root` | 29,148 | 0.01 |

log-Ne random-walk step scale tau = 5.181

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,198.0 | 222 | 49.39 |
| SNP | -36,026.9 | 6 | 12023.50 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +112.63 | -124.78 | -97.74 |
| **IBS** | -124.78 | +93.99 | +153.97 |
| **TSI** | -97.74 | +153.97 | +41.52 |

![spectrum](spectrum_fit.png)
