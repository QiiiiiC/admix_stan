# `(((IBS,TSI.1),EAS),TSI.2)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 19 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +40.62 | +- 1.15 (MC) |
| logZ (importance sampling) | +117.38 | |
| ESS of the IS weights | 1.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 7 | 24 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 135.4 +- 3.8 | 146.4 |
| 3 | MERGE | n1 + EAS -> n2 | 207.3 +- 2.9 | 353.6 |
| 4 | MERGE | TSI.1 + n2 -> root | 253.6 +- 1.5 | 607.2 |

## Admixture fraction

**f = 0.027 +- 0.003** (fraction from `TSI.1`; 0.973 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 21,264,532 | 16,045,736 |
| `IBS` | 3,388,340 | 2,351,966 |
| `TSI` | 2,250,483 | 1,419,834 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 880,281 | 0.04 |
| `IBS` | 276,559 | 0.05 |
| `TSI` | 468,503 | 0.12 |
| `TSI.1` | 989 | 0.27 |
| `TSI.2` | 549,817 | 0.16 |
| `n1` | 25,802 | 0.13 |
| `n2` | 39 | 0.06 |
| `root` | 11,538 | 0.00 |

log-Ne random-walk step scale tau_ibd = 2.415

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,535 | 1,625 |
| `IBS` | 342,910 | 352,349 |
| `TSI` | 57,932 | 58,124 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,060 | 0.02 |
| `IBS` | 397,758 | 0.03 |
| `TSI` | 56,661 | 0.03 |
| `TSI.1` | 8,732 | 0.05 |
| `TSI.2` | 62,676 | 0.02 |
| `n1` | 27,875 | 0.07 |
| `n2` | 28,674 | 0.09 |
| `root` | 18,631 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.668

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +333.3 | 222 | 36.31 |
| SNP | +23.8 | 6 | 6.60 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.06 | -0.61 | +0.73 |
| **IBS** | -0.61 | +1.28 | -0.13 |
| **TSI** | +0.73 | -0.13 | -1.28 |

![spectrum](spectrum_fit.png)
