# `(((EAS,IBS.1),IBS.2),TSI)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 10 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -3720.37 | +- 0.46 (MC) |
| logZ (importance sampling) | -3679.71 | |
| ESS of the IS weights | 1.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 6 | 24 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 88.9 +- 0.4 | 88.9 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 51.1 +- 0.5 | 140.0 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 32.0 +- 0.3 | 172.0 |
| 4 | MERGE | n2 + TSI -> root | 1.0 +- 0.0 | 173.0 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `IBS.1`; 0.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 5,131,979 | 4,499,342 |
| `IBS` | 1,283,207 | 1,068,848 |
| `TSI` | 1,647,782 | 1,311,290 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,383,191 | 0.07 |
| `IBS` | 781,794 | 0.07 |
| `TSI` | 306,110 | 0.03 |
| `IBS.1` | 38,767 | 0.03 |
| `IBS.2` | 16,499 | 0.03 |
| `n1` | 22,662 | 0.02 |
| `n2` | 122,666 | 0.01 |
| `root` | 90,965 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.850

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 745 | 710 |
| `IBS` | 172,612 | 169,099 |
| `TSI` | 100,340 | 106,292 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 807 | 0.02 |
| `IBS` | 170,679 | 0.03 |
| `TSI` | 110,288 | 0.05 |
| `IBS.1` | 124,203 | 0.03 |
| `IBS.2` | 4,607 | 0.02 |
| `n1` | 5,856 | 0.01 |
| `n2` | 17,866 | 0.01 |
| `root` | 18,416 | 0.01 |

log-Ne random-walk step scale tau_snp = 1.059

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -3,497.5 | 222 | 70.11 |
| SNP | +32.2 | 6 | 3.81 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.03 | -0.13 | +0.18 |
| **IBS** | -0.13 | +0.11 | +0.14 |
| **TSI** | +0.18 | +0.14 | -0.47 |

![spectrum](spectrum_fit.png)
