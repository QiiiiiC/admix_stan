# `(((EAS.1,TSI),IBS),EAS.2)`

**Normal, separate IBD/SNP Ne** | topology 07 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +3817.77 | +- 0.18 (MC) |
| logZ (importance sampling) | +3831.77 | |
| ESS of the IS weights | 9.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 1 / 11 | 22 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 125.6 +- 0.3 | 125.6 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 23.8 +- 0.7 | 149.4 |
| 3 | MERGE | n1 + IBS -> n2 | 1.5 +- 0.0 | 151.0 |
| 4 | MERGE | EAS.1 + n2 -> root | 136.7 +- 1.4 | 287.7 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,597,171 | 0.02 |
| `IBS` | 287,527 | 0.02 |
| `TSI` | 426,069 | 0.02 |
| `EAS.1` | 47,991 | 0.02 |
| `EAS.2` | 19,848 | 0.02 |
| `n1` | 20,524 | 0.02 |
| `n2` | 18,362 | 0.02 |
| `root` | 10,112 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.273

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,051 | 0.01 |
| `IBS` | 118,357 | 0.03 |
| `TSI` | 100,077 | 0.03 |
| `EAS.1` | 2,784 | 0.01 |
| `EAS.2` | 51,049 | 0.02 |
| `n1` | 50,971 | 0.02 |
| `n2` | 50,258 | 0.02 |
| `root` | 13,750 | 0.00 |

log-Ne random-walk step scale tau_snp = 0.814

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +3,961.4 | 222 | 3.63 |
| SNP | +41.5 | 6 | 0.69 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.03 | -0.03 | +0.09 |
| **IBS** | -0.03 | -0.04 | +0.10 |
| **TSI** | +0.09 | +0.10 | -0.26 |

![spectrum](spectrum_fit.png)
