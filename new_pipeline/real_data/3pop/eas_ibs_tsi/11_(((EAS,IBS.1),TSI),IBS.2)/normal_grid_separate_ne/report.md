# `(((EAS,IBS.1),TSI),IBS.2)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 11 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +483.43 | +- 0.13 (MC) |
| logZ (importance sampling) | +498.85 | |
| ESS of the IS weights | 1.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 7 | 24 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 127.8 +- 0.2 | 127.8 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 1.1 +- 0.0 | 128.8 |
| 3 | MERGE | n1 + TSI -> n2 | 15.6 +- 0.1 | 144.4 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.0 +- 0.0 | 145.4 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `IBS.1`; 0.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,345,487 | 3,845,379 |
| `IBS` | 2,579,764 | 2,011,770 |
| `TSI` | 1,933,118 | 1,511,933 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,569,505 | 0.03 |
| `IBS` | 575,749 | 0.03 |
| `TSI` | 316,824 | 0.03 |
| `IBS.1` | 6,626 | 0.03 |
| `IBS.2` | 25,055 | 0.02 |
| `n1` | 25,070 | 0.02 |
| `n2` | 170,311 | 0.01 |
| `root` | 81,673 | 0.01 |

log-Ne random-walk step scale tau_ibd = 2.034

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 799 | 789 |
| `IBS` | 188,986 | 178,669 |
| `TSI` | 140,953 | 136,579 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 716 | 0.01 |
| `IBS` | 129,739 | 0.03 |
| `TSI` | 93,435 | 0.04 |
| `IBS.1` | 28,085 | 0.02 |
| `IBS.2` | 9,822 | 0.01 |
| `n1` | 9,824 | 0.01 |
| `n2` | 32,235 | 0.01 |
| `root` | 31,037 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.800

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +745.9 | 222 | 31.74 |
| SNP | +39.0 | 6 | 1.54 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.10 | +0.02 | -0.21 |
| **IBS** | +0.02 | -0.71 | +0.72 |
| **TSI** | -0.21 | +0.72 | -0.27 |

![spectrum](spectrum_fit.png)
