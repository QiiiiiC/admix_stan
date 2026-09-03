# `(((IBS,TSI.1),TSI.2),EAS)`

**Normal, separate IBD/SNP Ne** | topology 20 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -97.13 | +- 0.81 (MC) |
| logZ (importance sampling) | -38.16 | |
| ESS of the IS weights | 1.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 7 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 13.5 +- 0.3 | 13.5 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 94.2 +- 0.6 | 107.7 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 51.4 +- 0.3 | 159.1 |
| 4 | MERGE | n2 + EAS -> root | 185.9 +- 0.6 | 345.0 |

## Admixture fraction

**f = 0.911 +- 0.001** (fraction from `TSI.1`; 0.089 from `TSI.2`)

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 966,633 | 0.04 |
| `IBS` | 637,521 | 0.06 |
| `TSI` | 1,364,355 | 0.04 |
| `TSI.1` | 292,017 | 0.04 |
| `TSI.2` | 339,080 | 0.05 |
| `n1` | 25,249 | 0.04 |
| `n2` | 32,829 | 0.02 |
| `root` | 49 | 0.03 |

log-Ne random-walk step scale tau_ibd = 1.681

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 105,702 | 0.06 |
| `IBS` | 146,159 | 0.08 |
| `TSI` | 174,639 | 0.08 |
| `TSI.1` | 121,699 | 0.07 |
| `TSI.2` | 107,154 | 0.07 |
| `n1` | 47,130 | 0.06 |
| `n2` | 1,050 | 0.02 |
| `root` | 3,474 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.992

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +135.4 | 222 | 37.76 |
| SNP | +21.8 | 6 | 7.27 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.13 | +0.55 | -0.29 |
| **IBS** | +0.55 | -0.95 | -0.10 |
| **TSI** | -0.29 | -0.10 | +0.65 |

![spectrum](spectrum_fit.png)
