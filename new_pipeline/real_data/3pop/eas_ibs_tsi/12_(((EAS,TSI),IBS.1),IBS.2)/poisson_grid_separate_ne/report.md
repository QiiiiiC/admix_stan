# `(((EAS,TSI),IBS.1),IBS.2)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 12 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -7025.04 | +- 0.58 (MC) |
| logZ (importance sampling) | -6978.40 | |
| ESS of the IS weights | 3.8 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 13 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 223.5 +- 3.1 | 223.5 |
| 2 | MERGE | EAS + TSI -> n1 | 1.0 +- 0.0 | 224.5 |
| 3 | MERGE | IBS.2 + n1 -> n2 | 1.0 +- 0.0 | 225.5 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.0 +- 0.0 | 226.5 |

## Admixture fraction

**f = 0.003 +- 0.002** (fraction from `IBS.1`; 0.997 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 87,340,835 | 49,072,397 |
| `IBS` | 5,551,454 | 3,594,364 |
| `TSI` | 3,678,308 | 2,325,152 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,053,577 | 0.03 |
| `IBS` | 403,044 | 0.07 |
| `TSI` | 322,268 | 0.04 |
| `IBS.1` | 9,906 | 0.13 |
| `IBS.2` | 50 | 0.11 |
| `n1` | 663 | 0.08 |
| `n2` | 1,262 | 0.10 |
| `root` | 9,598 | 0.12 |

log-Ne random-walk step scale tau_ibd = 3.379

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 2,440 | 2,163 |
| `IBS` | 245,702 | 226,842 |
| `TSI` | 269,054 | 247,751 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,222 | 0.03 |
| `IBS` | 162,664 | 0.14 |
| `TSI` | 177,447 | 0.14 |
| `IBS.1` | 19,195 | 0.18 |
| `IBS.2` | 11,456 | 0.58 |
| `n1` | 14,058 | 0.39 |
| `n2` | 15,115 | 0.34 |
| `root` | 19,015 | 0.18 |

log-Ne random-walk step scale tau_snp = 1.037

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -6,662.8 | 222 | 163.81 |
| SNP | +29.2 | 6 | 4.82 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.10 | +0.37 | -0.17 |
| **IBS** | +0.37 | -0.98 | +0.30 |
| **TSI** | -0.17 | +0.30 | +0.04 |

![spectrum](spectrum_fit.png)
