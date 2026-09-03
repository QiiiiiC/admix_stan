# `((EAS,TSI),IBS)`

**Poisson, shared Ne, recent grid** | topology 02 of 21 | tree (no admixture) | 2 events, 5 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -43444.45 | +- 0.28 (MC) |
| logZ (importance sampling) | -43427.99 | |
| ESS of the IS weights | 3.3 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -17.65 | already applied |
| seed kept / runtime | 1 | 3 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + TSI -> n1 | 301.4 +- 0.6 | 301.4 |
| 2 | MERGE | n1 + IBS -> root | 1.0 +- 0.0 | 302.4 |

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 661,728,743 | 187,834,651 |
| `IBS` | 1,432,764 | 993,220 |
| `TSI` | 623,230 | 538,419 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,050,945 | 0.02 |
| `IBS` | 232,787 | 0.05 |
| `TSI` | 315,452 | 0.03 |
| `n1` | 36 | 0.04 |
| `root` | 443 | 0.02 |

log-Ne random-walk step scale tau = 3.891

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,161.3 | 222 | 11738.31 |
| SNP | -36,030.3 | 6 | 12024.64 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +112.59 | -97.27 | -125.37 |
| **IBS** | -97.27 | +42.12 | +153.41 |
| **TSI** | -125.37 | +153.41 | +94.72 |

![spectrum](spectrum_fit.png)
