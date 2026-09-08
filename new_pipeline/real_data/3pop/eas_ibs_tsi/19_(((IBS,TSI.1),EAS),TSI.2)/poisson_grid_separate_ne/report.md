# `(((IBS,TSI.1),EAS),TSI.2)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 19 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -4421.94 | +- 0.13 (MC) |
| logZ (importance sampling) | -4407.74 | |
| ESS of the IS weights | 2.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 11 | 26 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 101.9 +- 0.5 | 101.9 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 9.5 +- 0.0 | 111.4 |
| 3 | MERGE | n1 + EAS -> n2 | 257.4 +- 0.6 | 368.8 |
| 4 | MERGE | TSI.1 + n2 -> root | 122.2 +- 0.7 | 491.1 |

## Admixture fraction

**f = 0.796 +- 0.002** (fraction from `TSI.1`; 0.204 from `TSI.2`)

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 57,743,856 | 33,451,048 |
| `IBS` | 4,733,574 | 3,148,680 |
| `TSI` | 1,281,531 | 878,967 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,002,452 | 0.02 |
| `IBS` | 657,912 | 0.04 |
| `TSI` | 409,645 | 0.03 |
| `TSI.1` | 60,331 | 0.03 |
| `TSI.2` | 72,416 | 0.04 |
| `n1` | 20,726 | 0.04 |
| `n2` | 7 | 0.03 |
| `root` | 10,232 | 0.00 |

log-Ne random-walk step scale tau_ibd = 2.509

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 2,506 | 2,306 |
| `IBS` | 434,419 | 410,583 |
| `TSI` | 649,019 | 575,726 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,041 | 0.01 |
| `IBS` | 376,453 | 0.03 |
| `TSI` | 646,756 | 0.04 |
| `TSI.1` | 433,111 | 0.05 |
| `TSI.2` | 284,436 | 0.03 |
| `n1` | 238,252 | 0.03 |
| `n2` | 58,671 | 0.00 |
| `root` | 58,713 | 0.00 |

log-Ne random-walk step scale tau_snp = 0.717

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,143.0 | 222 | 29.46 |
| SNP | +40.5 | 6 | 1.04 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.07 | -0.35 | +0.22 |
| **IBS** | -0.35 | +0.17 | +0.54 |
| **TSI** | +0.22 | +0.54 | -0.92 |

![spectrum](spectrum_fit.png)
