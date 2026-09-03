# `(((IBS.1,TSI),IBS.2),EAS)`

**Normal, separate IBD/SNP Ne** | topology 14 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -73.17 | +- 0.95 (MC) |
| logZ (importance sampling) | -13.24 | |
| ESS of the IS weights | 6.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 7 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 65.0 +- 0.3 | 65.0 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 61.8 +- 1.6 | 126.7 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 148.4 +- 10.6 | 275.1 |
| 4 | MERGE | n2 + EAS -> root | 62.1 +- 19.2 | 337.2 |

## Admixture fraction

**f = 0.072 +- 0.008** (fraction from `IBS.1`; 0.928 from `IBS.2`)

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 964,659 | 0.03 |
| `IBS` | 1,007,922 | 0.08 |
| `TSI` | 404,290 | 0.12 |
| `IBS.1` | 703 | 0.22 |
| `IBS.2` | 210,971 | 0.12 |
| `n1` | 72,213 | 0.12 |
| `n2` | 1,551 | 0.56 |
| `root` | 84 | 0.48 |

log-Ne random-walk step scale tau_ibd = 1.476

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 72,861 | 0.40 |
| `IBS` | 142,248 | 0.32 |
| `TSI` | 141,724 | 0.29 |
| `IBS.1` | 12,146 | 0.45 |
| `IBS.2` | 48,210 | 0.17 |
| `n1` | 22,778 | 0.15 |
| `n2` | 368 | 0.41 |
| `root` | 1,321 | 0.97 |

log-Ne random-walk step scale tau_snp = 1.098

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +131.3 | 222 | 37.92 |
| SNP | +33.3 | 6 | 3.43 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.18 | +0.18 | +0.18 |
| **IBS** | +0.18 | -0.24 | -0.11 |
| **TSI** | +0.18 | -0.11 | -0.23 |

![spectrum](spectrum_fit.png)
