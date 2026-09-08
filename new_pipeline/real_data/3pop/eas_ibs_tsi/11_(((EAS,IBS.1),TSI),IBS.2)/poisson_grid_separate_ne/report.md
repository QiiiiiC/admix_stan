# `(((EAS,IBS.1),TSI),IBS.2)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 11 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -3704.28 | +- 0.23 (MC) |
| logZ (importance sampling) | -3686.25 | |
| ESS of the IS weights | 5.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 4 | 26 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 87.7 +- 0.3 | 87.7 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 51.5 +- 0.4 | 139.2 |
| 3 | MERGE | n1 + TSI -> n2 | 34.8 +- 0.2 | 174.0 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.1 +- 0.0 | 175.1 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `IBS.1`; 0.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,219,881 | 3,911,237 |
| `IBS` | 1,774,990 | 1,520,927 |
| `TSI` | 2,128,671 | 1,435,537 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,407,347 | 0.05 |
| `IBS` | 757,955 | 0.05 |
| `TSI` | 305,177 | 0.02 |
| `IBS.1` | 40,744 | 0.02 |
| `IBS.2` | 33,564 | 0.06 |
| `n1` | 24,296 | 0.01 |
| `n2` | 145,514 | 0.01 |
| `root` | 83,886 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.905

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 745 | 792 |
| `IBS` | 162,390 | 162,527 |
| `TSI` | 120,223 | 114,483 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 794 | 0.01 |
| `IBS` | 172,357 | 0.04 |
| `TSI` | 117,417 | 0.05 |
| `IBS.1` | 98,184 | 0.04 |
| `IBS.2` | 8,951 | 0.03 |
| `n1` | 7,709 | 0.01 |
| `n2` | 19,581 | 0.01 |
| `root` | 20,700 | 0.01 |

log-Ne random-walk step scale tau_snp = 1.061

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -3,486.0 | 222 | 71.09 |
| SNP | +37.7 | 6 | 1.97 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.01 | +0.04 | -0.06 |
| **IBS** | +0.04 | -0.40 | +0.35 |
| **TSI** | -0.06 | +0.35 | -0.21 |

![spectrum](spectrum_fit.png)
