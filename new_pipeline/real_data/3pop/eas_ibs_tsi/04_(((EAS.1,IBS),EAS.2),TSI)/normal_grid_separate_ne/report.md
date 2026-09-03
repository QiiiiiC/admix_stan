# `(((EAS.1,IBS),EAS.2),TSI)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 04 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +446.84 | +- 0.40 (MC) |
| logZ (importance sampling) | +477.92 | |
| ESS of the IS weights | 2.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 7 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 134.6 +- 0.8 | 134.6 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 1.2 +- 0.0 | 135.8 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 7.6 +- 0.3 | 143.4 |
| 4 | MERGE | n2 + TSI -> root | 1.0 +- 0.0 | 144.4 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,639,957 | 3,999,167 |
| `IBS` | 2,885,045 | 2,166,558 |
| `TSI` | 1,802,394 | 1,276,223 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,541,070 | 0.05 |
| `IBS` | 540,012 | 0.07 |
| `TSI` | 314,756 | 0.02 |
| `EAS.1` | 13,020 | 0.04 |
| `EAS.2` | 2,706 | 0.06 |
| `n1` | 2,714 | 0.06 |
| `n2` | 137,063 | 0.01 |
| `root` | 85,651 | 0.01 |

log-Ne random-walk step scale tau_ibd = 2.027

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 758 | 753 |
| `IBS` | 147,441 | 144,839 |
| `TSI` | 113,435 | 111,251 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 748 | 0.02 |
| `IBS` | 135,060 | 0.04 |
| `TSI` | 101,721 | 0.04 |
| `EAS.1` | 15,795 | 0.01 |
| `EAS.2` | 22,358 | 0.01 |
| `n1` | 22,379 | 0.01 |
| `n2` | 26,084 | 0.01 |
| `root` | 25,397 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.669

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +725.4 | 222 | 31.90 |
| SNP | +35.3 | 6 | 2.76 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.02 | +0.03 | +0.01 |
| **IBS** | +0.03 | -0.15 | +0.10 |
| **TSI** | +0.01 | +0.10 | -0.11 |

![spectrum](spectrum_fit.png)
