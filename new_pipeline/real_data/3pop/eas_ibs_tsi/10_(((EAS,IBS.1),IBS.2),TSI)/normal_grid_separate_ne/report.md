# `(((EAS,IBS.1),IBS.2),TSI)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 10 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +498.83 | +- 0.13 (MC) |
| logZ (importance sampling) | +515.21 | |
| ESS of the IS weights | 1.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 7 | 27 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 92.1 +- 0.2 | 92.1 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 43.5 +- 0.2 | 135.6 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 9.1 +- 0.0 | 144.7 |
| 4 | MERGE | n2 + TSI -> root | 1.0 +- 0.0 | 145.7 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `IBS.1`; 0.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,540,289 | 4,220,429 |
| `IBS` | 1,788,885 | 1,211,100 |
| `TSI` | 1,845,464 | 1,209,027 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,510,428 | 0.04 |
| `IBS` | 759,086 | 0.04 |
| `TSI` | 316,907 | 0.03 |
| `IBS.1` | 32,831 | 0.02 |
| `IBS.2` | 10,074 | 0.02 |
| `n1` | 13,187 | 0.02 |
| `n2` | 131,191 | 0.01 |
| `root` | 80,896 | 0.01 |

log-Ne random-walk step scale tau_ibd = 2.039

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 663 | 681 |
| `IBS` | 167,894 | 165,938 |
| `TSI` | 110,281 | 110,313 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 762 | 0.01 |
| `IBS` | 165,089 | 0.03 |
| `TSI` | 113,885 | 0.02 |
| `IBS.1` | 87,735 | 0.02 |
| `IBS.2` | 11,899 | 0.01 |
| `n1` | 12,304 | 0.00 |
| `n2` | 17,999 | 0.00 |
| `root` | 18,724 | 0.00 |

log-Ne random-walk step scale tau_snp = 1.090

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +746.6 | 222 | 31.75 |
| SNP | +41.5 | 6 | 0.72 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.04 | +0.01 | +0.07 |
| **IBS** | +0.01 | +0.31 | -0.36 |
| **TSI** | +0.07 | -0.36 | +0.21 |

![spectrum](spectrum_fit.png)
