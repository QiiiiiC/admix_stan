# `((IBS,TSI),EAS)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 03 of 21 | tree (no admixture) | 2 events, 5 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5258.61 | +- 1.39 (MC) |
| logZ (importance sampling) | -5175.87 | |
| ESS of the IS weights | 1.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -26.08 | already applied |
| seed kept / runtime | 1 | 14 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | IBS + TSI -> n1 | 161.8 +- 2.5 | 161.8 |
| 2 | MERGE | EAS + n1 -> root | 167.1 +- 11.5 | 328.9 |

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 68,936,170 | 36,996,746 |
| `IBS` | 4,480,242 | 2,921,137 |
| `TSI` | 1,570,859 | 1,072,041 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 834,017 | 0.05 |
| `IBS` | 294,472 | 0.10 |
| `TSI` | 404,964 | 0.09 |
| `n1` | 13,676 | 0.07 |
| `root` | 634 | 0.45 |

log-Ne random-walk step scale tau_ibd = 2.548

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 54,344 | 45,906 |
| `IBS` | 218,681 | 197,395 |
| `TSI` | 133,946 | 122,521 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 21,341 | 0.19 |
| `IBS` | 117,802 | 0.11 |
| `TSI` | 97,080 | 0.12 |
| `n1` | 1,018 | 0.09 |
| `root` | 7,183 | 0.09 |

log-Ne random-walk step scale tau_snp = 1.305

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,012.0 | 222 | 41.37 |
| SNP | +31.9 | 6 | 3.90 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.01 | -0.09 | +0.08 |
| **IBS** | -0.09 | -0.24 | +0.45 |
| **TSI** | +0.08 | +0.45 | -0.58 |

![spectrum](spectrum_fit.png)
