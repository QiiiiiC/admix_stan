# `(((EAS,IBS),TSI.1),TSI.2)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 16 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -7324.66 | +- 0.23 (MC) |
| logZ (importance sampling) | -7305.36 | |
| ESS of the IS weights | 6.4 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 7 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 11.3 +- 0.0 | 11.3 |
| 2 | MERGE | EAS + IBS -> n1 | 228.4 +- 1.6 | 239.7 |
| 3 | MERGE | TSI.2 + n1 -> n2 | 1.0 +- 0.0 | 240.7 |
| 4 | MERGE | TSI.1 + n2 -> root | 1.0 +- 0.0 | 241.7 |

## Admixture fraction

**f = 0.006 +- 0.000** (fraction from `TSI.1`; 0.994 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 79,440,818 | 46,906,659 |
| `IBS` | 9,918,480 | 5,722,758 |
| `TSI` | 3,531,137 | 2,490,411 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,084,473 | 0.02 |
| `IBS` | 236,546 | 0.04 |
| `TSI` | 568,765 | 0.02 |
| `TSI.1` | 6,872 | 0.01 |
| `TSI.2` | 308,093 | 0.03 |
| `n1` | 286 | 0.02 |
| `n2` | 1,140 | 0.01 |
| `root` | 5,135 | 0.01 |

log-Ne random-walk step scale tau_ibd = 3.112

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,360 | 1,356 |
| `IBS` | 204,777 | 203,890 |
| `TSI` | 168,777 | 168,207 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,327 | 0.01 |
| `IBS` | 199,341 | 0.04 |
| `TSI` | 166,277 | 0.04 |
| `TSI.1` | 21,238 | 0.01 |
| `TSI.2` | 167,785 | 0.05 |
| `n1` | 20,884 | 0.01 |
| `n2` | 21,124 | 0.01 |
| `root` | 21,213 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.864

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,008.5 | 222 | 345.45 |
| SNP | +42.1 | 6 | 0.49 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.04 | -0.03 | +0.10 |
| **IBS** | -0.03 | +0.15 | -0.10 |
| **TSI** | +0.10 | -0.10 | -0.10 |

![spectrum](spectrum_fit.png)
