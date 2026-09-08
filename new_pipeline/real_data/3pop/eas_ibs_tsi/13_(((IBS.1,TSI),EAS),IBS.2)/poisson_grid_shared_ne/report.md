# `(((IBS.1,TSI),EAS),IBS.2)`

**Poisson, shared Ne, recent grid** | topology 13 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5094.04 | +- 0.44 (MC) |
| logZ (importance sampling) | -5062.77 | |
| ESS of the IS weights | 3.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 1 / 2 | 29 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 57.2 +- 1.4 | 57.2 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 194.0 +- 1.2 | 251.2 |
| 3 | MERGE | n1 + EAS -> n2 | 106.3 +- 3.8 | 357.5 |
| 4 | MERGE | IBS.1 + n2 -> root | 3.7 +- 0.0 | 361.1 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `IBS.1`; 1.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 60,493,828 | 31,442,340 |
| `IBS` | 1,488,649 | 1,369,311 |
| `TSI` | 1,532,295 | 1,084,672 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 829,497 | 0.02 |
| `IBS` | 1,027,407 | 0.14 |
| `TSI` | 389,409 | 0.05 |
| `IBS.1` | 193 | 0.09 |
| `IBS.2` | 115,421 | 0.09 |
| `n1` | 590 | 0.04 |
| `n2` | 304 | 0.18 |
| `root` | 181 | 0.08 |

log-Ne random-walk step scale tau = 2.252

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,864.6 | 222 | 156.02 |
| SNP | +32.8 | 6 | 3.61 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.05 | +0.55 | -0.65 |
| **IBS** | +0.55 | -0.63 | -0.44 |
| **TSI** | -0.65 | -0.44 | +1.66 |

![spectrum](spectrum_fit.png)
