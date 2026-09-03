# `((EAS,TSI),IBS)`

**Normal, shared Ne, recent grid** | topology 02 of 21 | tree (no admixture) | 2 events, 5 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -40979.51 | +- 0.39 (MC) |
| logZ (importance sampling) | -40956.08 | |
| ESS of the IS weights | 2.3 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -17.65 | already applied |
| seed kept / runtime | 7 | 8 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + TSI -> n1 | 162.9 +- 0.6 | 162.9 |
| 2 | MERGE | n1 + IBS -> root | 2.8 +- 0.2 | 165.7 |

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 55,293,808 | 27,760,069 |
| `IBS` | 3,562,701 | 2,268,685 |
| `TSI` | 1,389,556 | 1,068,282 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,412,884 | 0.03 |
| `IBS` | 211,334 | 0.04 |
| `TSI` | 334,319 | 0.05 |
| `n1` | 7,747 | 0.11 |
| `root` | 29,282 | 0.03 |

log-Ne random-walk step scale tau = 2.444

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,235.0 | 222 | 49.65 |
| SNP | -39,604.6 | 6 | 13216.07 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +117.34 | -115.18 | -116.74 |
| **IBS** | -115.18 | +110.89 | +116.26 |
| **TSI** | -116.74 | +116.26 | +113.22 |

![spectrum](spectrum_fit.png)
