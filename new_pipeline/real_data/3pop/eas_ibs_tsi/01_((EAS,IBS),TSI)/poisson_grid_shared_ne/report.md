# `((EAS,IBS),TSI)`

**Poisson, shared Ne, recent grid** | topology 01 of 21 | tree (no admixture) | 2 events, 5 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -43210.03 | +- 0.99 (MC) |
| logZ (importance sampling) | -43168.61 | |
| ESS of the IS weights | 11.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -17.65 | already applied |
| mode kept / MAP start / runtime | 2 / 6 | 17 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + IBS -> n1 | 294.1 +- 1.8 | 294.1 |
| 2 | MERGE | n1 + TSI -> root | 1.4 +- 0.0 | 295.5 |

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 53,100,324 | 34,036,860 |
| `IBS` | 2,906,935 | 2,189,453 |
| `TSI` | 2,168,164 | 1,511,847 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,076,503 | 0.03 |
| `IBS` | 229,518 | 0.06 |
| `TSI` | 309,859 | 0.07 |
| `n1` | 50 | 0.05 |
| `root` | 614 | 0.09 |

log-Ne random-walk step scale tau = 3.088

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,096.9 | 222 | 8381.47 |
| SNP | -35,890.1 | 6 | 11977.92 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +112.48 | -124.32 | -97.91 |
| **IBS** | -124.32 | +92.70 | +154.40 |
| **TSI** | -97.91 | +154.40 | +41.44 |

![spectrum](spectrum_fit.png)
