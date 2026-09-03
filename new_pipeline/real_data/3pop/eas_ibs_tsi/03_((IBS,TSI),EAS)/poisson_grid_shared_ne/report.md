# `((IBS,TSI),EAS)`

**Poisson, shared Ne, recent grid** | topology 03 of 21 | tree (no admixture) | 2 events, 5 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5249.68 | +- 0.41 (MC) |
| logZ (importance sampling) | -5223.03 | |
| ESS of the IS weights | 3.8 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -17.65 | already applied |
| seed kept / runtime | 13 | 12 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | IBS + TSI -> n1 | 247.7 +- 1.0 | 247.7 |
| 2 | MERGE | EAS + n1 -> root | 108.1 +- 1.5 | 355.8 |

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 80,616,327 | 38,695,291 |
| `IBS` | 4,480,425 | 2,861,320 |
| `TSI` | 1,609,473 | 1,094,528 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 829,934 | 0.03 |
| `IBS` | 293,115 | 0.03 |
| `TSI` | 398,372 | 0.05 |
| `n1` | 599 | 0.02 |
| `root` | 202 | 0.04 |

log-Ne random-walk step scale tau = 2.493

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,082.8 | 222 | 125.99 |
| SNP | +30.6 | 6 | 4.33 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.21 | +0.37 | +0.05 |
| **IBS** | +0.37 | +1.27 | -2.11 |
| **TSI** | +0.05 | -2.11 | +1.90 |

![spectrum](spectrum_fit.png)
