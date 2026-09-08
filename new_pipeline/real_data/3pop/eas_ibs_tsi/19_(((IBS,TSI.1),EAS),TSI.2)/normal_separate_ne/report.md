# `(((IBS,TSI.1),EAS),TSI.2)`

**Normal, separate IBD/SNP Ne** | topology 19 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +924.64 | +- 0.32 (MC) |
| logZ (importance sampling) | +950.02 | |
| ESS of the IS weights | 6.8 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 3 / 0 | 26 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 16.5 +- 0.1 | 16.5 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 100.2 +- 0.7 | 116.6 |
| 3 | MERGE | n1 + EAS -> n2 | 214.7 +- 1.8 | 331.4 |
| 4 | MERGE | TSI.1 + n2 -> root | 191.1 +- 4.3 | 522.5 |

## Admixture fraction

**f = 0.805 +- 0.003** (fraction from `TSI.1`; 0.195 from `TSI.2`)

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,205,385 | 0.03 |
| `IBS` | 691,646 | 0.04 |
| `TSI` | 2,375,651 | 0.05 |
| `TSI.1` | 185,379 | 0.05 |
| `TSI.2` | 2,386,931 | 0.07 |
| `n1` | 17,671 | 0.03 |
| `n2` | 44 | 0.08 |
| `root` | 11,174 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.716

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,840 | 0.01 |
| `IBS` | 261,016 | 0.03 |
| `TSI` | 432,062 | 0.03 |
| `TSI.1` | 470,266 | 0.03 |
| `TSI.2` | 255,313 | 0.03 |
| `n1` | 203,255 | 0.03 |
| `n2` | 125,477 | 0.01 |
| `root` | 25,626 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.865

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +1,098.6 | 222 | 29.18 |
| SNP | +36.6 | 6 | 2.36 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.04 | -0.11 | +0.03 |
| **IBS** | -0.11 | -0.05 | +0.27 |
| **TSI** | +0.03 | +0.27 | -0.32 |

![spectrum](spectrum_fit.png)
