# `((EAS,IBS.1),(IBS.2,TSI))`

**Poisson, separate IBD/SNP Ne** | topology 15 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -523.59 | +- 0.36 (MC) |
| logZ (importance sampling) | -495.91 | |
| ESS of the IS weights | 1.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 1 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 76.3 +- 0.6 | 76.3 |
| 2 | MERGE | IBS.1 + EAS -> n1 | 58.2 +- 0.5 | 134.5 |
| 3 | MERGE | IBS.2 + TSI -> n2 | 1.3 +- 0.0 | 135.8 |
| 4 | MERGE | n1 + n2 -> root | 131.0 +- 0.7 | 266.8 |

## Admixture fraction

**f = 0.006 +- 0.000** (fraction from `IBS.1`; 0.994 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,491,628 | 0.05 |
| `IBS` | 883,598 | 0.06 |
| `TSI` | 408,771 | 0.03 |
| `IBS.1` | 36,141 | 0.03 |
| `IBS.2` | 61,616 | 0.04 |
| `n1` | 36,063 | 0.03 |
| `n2` | 34,850 | 0.02 |
| `root` | 28,279 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.307

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,156 | 0.01 |
| `IBS` | 61,155 | 0.05 |
| `TSI` | 143,291 | 0.04 |
| `IBS.1` | 2,123 | 0.01 |
| `IBS.2` | 77,540 | 0.04 |
| `n1` | 2,126 | 0.01 |
| `n2` | 65,494 | 0.03 |
| `root` | 29,014 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.656

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -376.8 | 222 | 1.05 |
| SNP | +36.1 | 6 | 2.50 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.08 | +0.10 | -0.27 |
| **IBS** | +0.10 | -0.61 | +0.45 |
| **TSI** | -0.27 | +0.45 | +0.08 |

![spectrum](spectrum_fit.png)
