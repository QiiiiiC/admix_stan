# `(((IBS.1,TSI),EAS),IBS.2)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 13 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5267.52 | +- 0.48 (MC) |
| logZ (importance sampling) | -5232.14 | |
| ESS of the IS weights | 6.4 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 1 | 24 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 156.3 +- 0.5 | 156.3 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 1.0 +- 0.0 | 157.3 |
| 3 | MERGE | n1 + EAS -> n2 | 167.1 +- 2.9 | 324.5 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.1 +- 0.0 | 325.6 |

## Admixture fraction

**f = 0.002 +- 0.000** (fraction from `IBS.1`; 0.998 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 69,139,594 | 37,691,097 |
| `IBS` | 4,528,111 | 2,905,263 |
| `TSI` | 1,570,172 | 1,067,458 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 830,390 | 0.03 |
| `IBS` | 296,944 | 0.04 |
| `TSI` | 401,799 | 0.06 |
| `IBS.1` | 561 | 0.01 |
| `IBS.2` | 10,869 | 0.03 |
| `n1` | 15,885 | 0.03 |
| `n2` | 741 | 0.01 |
| `root` | 670 | 0.01 |

log-Ne random-walk step scale tau_ibd = 2.544

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 3,530 | 3,229 |
| `IBS` | 151,508 | 144,535 |
| `TSI` | 126,928 | 124,235 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,184 | 0.01 |
| `IBS` | 119,540 | 0.04 |
| `TSI` | 117,702 | 0.03 |
| `IBS.1` | 10,665 | 0.00 |
| `IBS.2` | 5,157 | 0.02 |
| `n1` | 4,977 | 0.02 |
| `n2` | 10,918 | 0.00 |
| `root` | 10,832 | 0.00 |

log-Ne random-walk step scale tau_snp = 0.838

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,954.9 | 222 | 40.67 |
| SNP | +31.6 | 6 | 4.00 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.10 | +0.15 | +0.04 |
| **IBS** | +0.15 | -0.35 | +0.06 |
| **TSI** | +0.04 | +0.06 | -0.13 |

![spectrum](spectrum_fit.png)
