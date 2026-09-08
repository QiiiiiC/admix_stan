# `(((EAS,TSI.1),IBS),TSI.2)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 17 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -4231.97 | +- 0.82 (MC) |
| logZ (importance sampling) | -4169.86 | |
| ESS of the IS weights | 1.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 7 | 26 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 132.8 +- 2.1 | 132.8 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 3.5 +- 0.4 | 136.3 |
| 3 | MERGE | n1 + IBS -> n2 | 43.7 +- 2.5 | 180.0 |
| 4 | MERGE | TSI.1 + n2 -> root | 1.0 +- 0.0 | 181.1 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,529,178 | 3,993,407 |
| `IBS` | 4,763,690 | 2,559,860 |
| `TSI` | 1,212,126 | 899,675 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,403,915 | 0.06 |
| `IBS` | 223,293 | 0.04 |
| `TSI` | 396,952 | 0.08 |
| `TSI.1` | 33,864 | 0.10 |
| `TSI.2` | 31,657 | 0.07 |
| `n1` | 29,910 | 0.09 |
| `n2` | 100,869 | 0.07 |
| `root` | 67,222 | 0.03 |

log-Ne random-walk step scale tau_ibd = 1.924

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 738 | 747 |
| `IBS` | 132,252 | 136,618 |
| `TSI` | 178,988 | 175,297 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 785 | 0.02 |
| `IBS` | 134,442 | 0.16 |
| `TSI` | 165,424 | 0.18 |
| `TSI.1` | 80,439 | 0.15 |
| `TSI.2` | 7,560 | 0.03 |
| `n1` | 7,194 | 0.04 |
| `n2` | 22,762 | 0.05 |
| `root` | 23,136 | 0.06 |

log-Ne random-walk step scale tau_snp = 1.017

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -3,982.2 | 222 | 80.06 |
| SNP | +27.5 | 6 | 5.38 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.05 | +0.06 | +0.04 |
| **IBS** | +0.06 | -0.27 | +0.17 |
| **TSI** | +0.04 | +0.17 | -0.24 |

![spectrum](spectrum_fit.png)
