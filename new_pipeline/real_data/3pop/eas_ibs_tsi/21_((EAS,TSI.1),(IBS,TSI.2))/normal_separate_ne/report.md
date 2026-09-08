# `((EAS,TSI.1),(IBS,TSI.2))`

**Normal, separate IBD/SNP Ne** | topology 21 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +3808.05 | +- 0.32 (MC) |
| logZ (importance sampling) | +3835.26 | |
| ESS of the IS weights | 1.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 1 / 2 | 22 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 39.2 +- 2.7 | 39.2 |
| 2 | MERGE | TSI.1 + EAS -> n1 | 89.0 +- 3.2 | 128.1 |
| 3 | MERGE | TSI.2 + IBS -> n2 | 22.9 +- 1.9 | 151.0 |
| 4 | MERGE | n1 + n2 -> root | 130.6 +- 2.1 | 281.6 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `TSI.1`; 1.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,560,445 | 0.04 |
| `IBS` | 289,779 | 0.04 |
| `TSI` | 405,925 | 0.05 |
| `TSI.1` | 74,246 | 0.08 |
| `TSI.2` | 444,194 | 0.04 |
| `n1` | 44,412 | 0.04 |
| `n2` | 18,036 | 0.05 |
| `root` | 12,728 | 0.10 |

log-Ne random-walk step scale tau_ibd = 1.273

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,028 | 0.01 |
| `IBS` | 115,062 | 0.06 |
| `TSI` | 107,762 | 0.06 |
| `TSI.1` | 2,611 | 0.01 |
| `TSI.2` | 99,558 | 0.06 |
| `n1` | 2,886 | 0.02 |
| `n2` | 50,471 | 0.03 |
| `root` | 14,267 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.834

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +3,948.5 | 222 | 3.64 |
| SNP | +39.7 | 6 | 1.31 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.00 | -0.01 | +0.02 |
| **IBS** | -0.01 | -0.10 | +0.14 |
| **TSI** | +0.02 | +0.14 | -0.18 |

![spectrum](spectrum_fit.png)
