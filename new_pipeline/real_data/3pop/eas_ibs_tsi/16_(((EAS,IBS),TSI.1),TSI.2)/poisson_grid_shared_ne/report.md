# `(((EAS,IBS),TSI.1),TSI.2)`

**Poisson, shared Ne, recent grid** | topology 16 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -43238.53 | +- 0.34 (MC) |
| logZ (importance sampling) | -43215.57 | |
| ESS of the IS weights | 3.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 13 | 22 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | EAS + IBS -> n1 | 280.8 +- 0.8 | 291.8 |
| 3 | MERGE | TSI.2 + n1 -> n2 | 1.0 +- 0.0 | 292.8 |
| 4 | MERGE | TSI.1 + n2 -> root | 1.0 +- 0.0 | 293.8 |

## Admixture fraction

**f = 0.001 +- 0.000** (fraction from `TSI.1`; 0.999 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 49,532,070 | 31,036,240 |
| `IBS` | 7,639,869 | 4,561,712 |
| `TSI` | 2,723,954 | 1,872,890 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,085,829 | 0.03 |
| `IBS` | 228,342 | 0.02 |
| `TSI` | 449,736 | 0.04 |
| `TSI.1` | 917 | 0.08 |
| `TSI.2` | 312,161 | 0.03 |
| `n1` | 35 | 0.04 |
| `n2` | 147 | 0.03 |
| `root` | 814 | 0.01 |

log-Ne random-walk step scale tau = 2.754

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,066.1 | 222 | 6069.11 |
| SNP | -35,887.5 | 6 | 11977.03 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +112.54 | -124.22 | -98.11 |
| **IBS** | -124.22 | +92.92 | +153.97 |
| **TSI** | -98.11 | +153.97 | +42.23 |

![spectrum](spectrum_fit.png)
