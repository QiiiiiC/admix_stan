# `(((EAS.1,IBS),TSI),EAS.2)`

**Poisson, separate IBD/SNP Ne** | topology 05 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -529.32 | +- 0.29 (MC) |
| logZ (importance sampling) | -506.61 | |
| ESS of the IS weights | 2.4 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 2 / 6 | 19 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 122.7 +- 0.7 | 122.7 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 1.9 +- 0.0 | 124.6 |
| 3 | MERGE | n1 + TSI -> n2 | 20.5 +- 0.3 | 145.2 |
| 4 | MERGE | EAS.1 + n2 -> root | 143.7 +- 2.7 | 288.8 |

## Admixture fraction

**f = 0.999 +- 0.000** (fraction from `EAS.1`; 0.001 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,642,791 | 0.03 |
| `IBS` | 623,241 | 0.03 |
| `TSI` | 408,127 | 0.06 |
| `EAS.1` | 52,014 | 0.04 |
| `EAS.2` | 10,297 | 0.02 |
| `n1` | 10,142 | 0.02 |
| `n2` | 26,099 | 0.02 |
| `root` | 12,925 | 0.11 |

log-Ne random-walk step scale tau_ibd = 1.362

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,197 | 0.02 |
| `IBS` | 119,915 | 0.12 |
| `TSI` | 99,970 | 0.08 |
| `EAS.1` | 2,199 | 0.02 |
| `EAS.2` | 61,975 | 0.06 |
| `n1` | 62,717 | 0.06 |
| `n2` | 56,677 | 0.05 |
| `root` | 13,537 | 0.03 |

log-Ne random-walk step scale tau_snp = 0.824

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -386.2 | 222 | 1.18 |
| SNP | +37.2 | 6 | 2.13 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.00 | +0.02 | -0.03 |
| **IBS** | +0.02 | -0.22 | +0.20 |
| **TSI** | -0.03 | +0.20 | -0.13 |

![spectrum](spectrum_fit.png)
