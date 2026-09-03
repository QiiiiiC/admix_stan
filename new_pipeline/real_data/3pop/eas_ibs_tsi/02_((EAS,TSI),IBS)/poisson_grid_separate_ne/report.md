# `((EAS,TSI),IBS)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 02 of 21 | tree (no admixture) | 2 events, 5 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -7332.78 | +- 0.66 (MC) |
| logZ (importance sampling) | -7296.48 | |
| ESS of the IS weights | 5.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -26.08 | already applied |
| seed kept / runtime | 1 | 14 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + TSI -> n1 | 233.4 +- 0.7 | 233.4 |
| 2 | MERGE | n1 + IBS -> root | 1.0 +- 0.0 | 234.4 |

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 67,148,022 | 34,401,212 |
| `IBS` | 11,512,324 | 6,105,589 |
| `TSI` | 5,231,049 | 2,939,477 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,029,048 | 0.01 |
| `IBS` | 231,016 | 0.06 |
| `TSI` | 320,887 | 0.10 |
| `n1` | 970 | 0.01 |
| `root` | 5,148 | 0.01 |

log-Ne random-walk step scale tau_ibd = 2.888

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 3,123 | 2,616 |
| `IBS` | 235,267 | 223,353 |
| `TSI` | 181,672 | 177,874 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,264 | 0.02 |
| `IBS` | 181,918 | 0.17 |
| `TSI` | 164,793 | 0.15 |
| `n1` | 9,704 | 0.02 |
| `root` | 13,242 | 0.02 |

log-Ne random-walk step scale tau_snp = 0.973

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,082.9 | 222 | 265.19 |
| SNP | +33.1 | 6 | 3.49 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.03 | +0.02 | -0.08 |
| **IBS** | +0.02 | -0.27 | +0.25 |
| **TSI** | -0.08 | +0.25 | -0.08 |

![spectrum](spectrum_fit.png)
