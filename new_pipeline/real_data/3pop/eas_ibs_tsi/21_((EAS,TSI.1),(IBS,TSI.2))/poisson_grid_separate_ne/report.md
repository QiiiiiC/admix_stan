# `((EAS,TSI.1),(IBS,TSI.2))`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 21 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -812.57 | +- 0.22 (MC) |
| logZ (importance sampling) | -792.32 | |
| ESS of the IS weights | 1.3 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 2 / 8 | 25 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 17.4 +- 0.1 | 17.4 |
| 2 | MERGE | TSI.1 + EAS -> n1 | 113.3 +- 0.3 | 130.7 |
| 3 | MERGE | TSI.2 + IBS -> n2 | 33.1 +- 0.3 | 163.8 |
| 4 | MERGE | n1 + n2 -> root | 110.2 +- 1.0 | 274.0 |

## Admixture fraction

**f = 0.001 +- 0.000** (fraction from `TSI.1`; 0.999 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,119,046 | 4,092,270 |
| `IBS` | 3,124,874 | 1,901,807 |
| `TSI` | 1,089,498 | 806,047 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,405,702 | 0.03 |
| `IBS` | 294,810 | 0.03 |
| `TSI` | 435,709 | 0.05 |
| `TSI.1` | 25,873 | 0.04 |
| `TSI.2` | 403,998 | 0.04 |
| `n1` | 41,342 | 0.01 |
| `n2` | 11,976 | 0.01 |
| `root` | 17,599 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.845

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,152 | 1,162 |
| `IBS` | 136,841 | 131,386 |
| `TSI` | 116,400 | 116,738 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,104 | 0.02 |
| `IBS` | 123,536 | 0.04 |
| `TSI` | 104,113 | 0.04 |
| `TSI.1` | 2,375 | 0.01 |
| `TSI.2` | 107,326 | 0.04 |
| `n1` | 2,368 | 0.01 |
| `n2` | 57,569 | 0.03 |
| `root` | 11,822 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.846

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -600.2 | 222 | 3.41 |
| SNP | +35.6 | 6 | 2.67 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.00 | +0.01 | -0.01 |
| **IBS** | +0.01 | -0.19 | +0.19 |
| **TSI** | -0.01 | +0.19 | -0.16 |

![spectrum](spectrum_fit.png)
