# `((EAS,IBS),TSI)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 01 of 21 | tree (no admixture) | 2 events, 5 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1360.58 | +- 0.15 (MC) |
| logZ (importance sampling) | -1344.01 | |
| ESS of the IS weights | 1.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -26.08 | already applied |
| mode kept / MAP start / runtime | 1 / 0 | 17 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + IBS -> n1 | 157.7 +- 0.3 | 157.7 |
| 2 | MERGE | n1 + TSI -> root | 4.1 +- 0.1 | 161.8 |

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 27,206,182 | 16,573,135 |
| `IBS` | 2,775,643 | 1,779,730 |
| `TSI` | 1,538,975 | 1,097,994 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,389,526 | 0.02 |
| `IBS` | 232,343 | 0.02 |
| `TSI` | 337,483 | 0.03 |
| `n1` | 17,125 | 0.03 |
| `root` | 32,143 | 0.03 |

log-Ne random-walk step scale tau_ibd = 2.417

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 860 | 865 |
| `IBS` | 111,501 | 111,093 |
| `TSI` | 125,895 | 125,497 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 876 | 0.01 |
| `IBS` | 110,278 | 0.03 |
| `TSI` | 124,710 | 0.03 |
| `n1` | 18,088 | 0.01 |
| `root` | 19,585 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.899

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,176.2 | 222 | 49.07 |
| SNP | +39.5 | 6 | 1.36 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.14 | -0.12 | -0.15 |
| **IBS** | -0.12 | -0.23 | +0.49 |
| **TSI** | -0.15 | +0.49 | -0.17 |

![spectrum](spectrum_fit.png)
