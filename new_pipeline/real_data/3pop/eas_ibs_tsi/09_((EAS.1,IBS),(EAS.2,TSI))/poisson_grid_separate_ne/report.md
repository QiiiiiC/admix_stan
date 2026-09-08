# `((EAS.1,IBS),(EAS.2,TSI))`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 09 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -3760.59 | +- 0.38 (MC) |
| logZ (importance sampling) | -3733.30 | |
| ESS of the IS weights | 1.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 2 | 25 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 133.0 +- 0.7 | 133.0 |
| 2 | MERGE | EAS.1 + IBS -> n1 | 2.1 +- 0.0 | 135.1 |
| 3 | MERGE | EAS.2 + TSI -> n2 | 34.7 +- 0.5 | 169.8 |
| 4 | MERGE | n1 + n2 -> root | 1.0 +- 0.0 | 170.9 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `EAS.1`; 1.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,138,508 | 3,586,772 |
| `IBS` | 1,696,141 | 1,368,070 |
| `TSI` | 1,851,448 | 1,293,140 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,492,997 | 0.04 |
| `IBS` | 550,528 | 0.07 |
| `TSI` | 306,074 | 0.06 |
| `EAS.1` | 7,536 | 0.05 |
| `EAS.2` | 30,200 | 0.03 |
| `n1` | 8,010 | 0.05 |
| `n2` | 149,621 | 0.02 |
| `root` | 97,434 | 0.02 |

log-Ne random-walk step scale tau_ibd = 1.757

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 748 | 762 |
| `IBS` | 194,653 | 187,577 |
| `TSI` | 95,512 | 104,326 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 759 | 0.01 |
| `IBS` | 188,312 | 0.03 |
| `TSI` | 107,517 | 0.05 |
| `EAS.1` | 57,230 | 0.02 |
| `EAS.2` | 7,533 | 0.00 |
| `n1` | 55,022 | 0.02 |
| `n2` | 19,695 | 0.00 |
| `root` | 19,461 | 0.00 |

log-Ne random-walk step scale tau_snp = 1.050

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -3,513.9 | 222 | 68.41 |
| SNP | +39.4 | 6 | 1.42 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.06 | -0.08 | -0.04 |
| **IBS** | -0.08 | -0.31 | +0.50 |
| **TSI** | -0.04 | +0.50 | -0.39 |

![spectrum](spectrum_fit.png)
