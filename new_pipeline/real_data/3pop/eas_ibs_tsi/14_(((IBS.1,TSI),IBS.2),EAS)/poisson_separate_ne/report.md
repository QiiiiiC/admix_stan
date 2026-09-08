# `(((IBS.1,TSI),IBS.2),EAS)`

**Poisson, separate IBD/SNP Ne** | topology 14 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5515.21 | +- 0.51 (MC) |
| logZ (importance sampling) | -5483.48 | |
| ESS of the IS weights | 3.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 3 / 7 | 18 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 72.9 +- 0.7 | 72.9 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 56.8 +- 0.8 | 129.7 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 15.9 +- 0.2 | 145.6 |
| 4 | MERGE | n2 + EAS -> root | 178.9 +- 1.8 | 324.6 |

## Admixture fraction

**f = 0.998 +- 0.000** (fraction from `IBS.1`; 0.002 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 898,711 | 0.04 |
| `IBS` | 910,586 | 0.07 |
| `TSI` | 422,981 | 0.07 |
| `IBS.1` | 71,247 | 0.04 |
| `IBS.2` | 66,714 | 0.04 |
| `n1` | 93,008 | 0.04 |
| `n2` | 27,804 | 0.03 |
| `root` | 511 | 0.02 |

log-Ne random-walk step scale tau_ibd = 1.250

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,815 | 0.01 |
| `IBS` | 99,754 | 0.06 |
| `TSI` | 104,159 | 0.06 |
| `IBS.1` | 133,079 | 0.06 |
| `IBS.2` | 137,061 | 0.05 |
| `n1` | 124,214 | 0.05 |
| `n2` | 113,591 | 0.05 |
| `root` | 18,467 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.679

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,364.9 | 222 | 43.64 |
| SNP | +36.1 | 6 | 2.52 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.05 | +0.06 | +0.03 |
| **IBS** | +0.06 | +0.01 | -0.13 |
| **TSI** | +0.03 | -0.13 | +0.06 |

![spectrum](spectrum_fit.png)
