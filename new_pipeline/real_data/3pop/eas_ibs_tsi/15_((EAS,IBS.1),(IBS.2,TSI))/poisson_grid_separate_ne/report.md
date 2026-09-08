# `((EAS,IBS.1),(IBS.2,TSI))`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 15 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -549.43 | +- 0.44 (MC) |
| logZ (importance sampling) | -518.55 | |
| ESS of the IS weights | 2.9 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 3 | 25 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 72.1 +- 1.6 | 72.1 |
| 2 | MERGE | IBS.1 + EAS -> n1 | 58.7 +- 1.8 | 130.8 |
| 3 | MERGE | IBS.2 + TSI -> n2 | 20.3 +- 0.5 | 151.1 |
| 4 | MERGE | n1 + n2 -> root | 139.5 +- 0.5 | 290.6 |

## Admixture fraction

**f = 0.006 +- 0.000** (fraction from `IBS.1`; 0.994 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,827,530 | 3,883,861 |
| `IBS` | 1,371,458 | 1,298,804 |
| `TSI` | 1,027,199 | 807,863 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,413,233 | 0.05 |
| `IBS` | 860,905 | 0.06 |
| `TSI` | 393,768 | 0.02 |
| `IBS.1` | 28,623 | 0.06 |
| `IBS.2` | 74,709 | 0.05 |
| `n1` | 41,326 | 0.04 |
| `n2` | 21,259 | 0.02 |
| `root` | 11,530 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.506

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,190 | 1,095 |
| `IBS` | 89,132 | 89,786 |
| `TSI` | 127,689 | 135,451 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,204 | 0.01 |
| `IBS` | 88,158 | 0.03 |
| `TSI` | 135,533 | 0.03 |
| `IBS.1` | 2,453 | 0.01 |
| `IBS.2` | 84,552 | 0.03 |
| `n1` | 2,312 | 0.01 |
| `n2` | 57,626 | 0.02 |
| `root` | 15,678 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.828

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -349.8 | 222 | 1.09 |
| SNP | +37.2 | 6 | 2.13 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.06 | -0.09 | -0.04 |
| **IBS** | -0.09 | -0.10 | +0.28 |
| **TSI** | -0.04 | +0.28 | -0.19 |

![spectrum](spectrum_fit.png)
