# `(((IBS,TSI.1),TSI.2),EAS)`

**Normal, shared Ne** | topology 20 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -240.14 | +- 1.57 (MC) |
| logZ (importance sampling) | -121.43 | |
| ESS of the IS weights | 1.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 7 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 18.8 +- 1.2 | 18.8 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 83.8 +- 3.5 | 102.6 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 253.3 +- 2.6 | 355.9 |
| 4 | MERGE | n2 + EAS -> root | 9.8 +- 2.6 | 365.7 |

## Admixture fraction

**f = 0.752 +- 0.013** (fraction from `TSI.1`; 0.248 from `TSI.2`)

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 977,262 | 0.02 |
| `IBS` | 609,976 | 0.18 |
| `TSI` | 916,503 | 0.10 |
| `TSI.1` | 185,563 | 0.07 |
| `TSI.2` | 421,503 | 0.20 |
| `n1` | 35,915 | 0.15 |
| `n2` | 55 | 0.36 |
| `root` | 7 | 0.32 |

log-Ne random-walk step scale tau = 1.041

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +54.1 | 222 | 38.34 |
| SNP | -23.0 | 6 | 22.22 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.65 | +2.24 | -3.54 |
| **IBS** | +2.24 | -8.43 | +4.47 |
| **TSI** | -3.54 | +4.47 | +2.55 |

![spectrum](spectrum_fit.png)
