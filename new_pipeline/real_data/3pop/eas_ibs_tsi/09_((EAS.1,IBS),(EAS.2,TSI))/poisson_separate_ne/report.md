# `((EAS.1,IBS),(EAS.2,TSI))`

**Poisson, separate IBD/SNP Ne** | topology 09 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -7174.96 | +- 0.13 (MC) |
| logZ (importance sampling) | -7162.31 | |
| ESS of the IS weights | 5.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 1 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | EAS.1 + IBS -> n1 | 79.4 +- 0.2 | 80.4 |
| 3 | MERGE | EAS.2 + TSI -> n2 | 149.1 +- 0.4 | 229.5 |
| 4 | MERGE | n1 + n2 -> root | 8.0 +- 0.0 | 237.5 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `EAS.1`; 1.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,199,571 | 0.02 |
| `IBS` | 863,223 | 0.09 |
| `TSI` | 344,113 | 0.02 |
| `EAS.1` | 46,076 | 0.02 |
| `EAS.2` | 1,198,579 | 0.02 |
| `n1` | 52,445 | 0.02 |
| `n2` | 2,299 | 0.02 |
| `root` | 5,095 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.439

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,260 | 0.01 |
| `IBS` | 538,197 | 0.04 |
| `TSI` | 145,811 | 0.03 |
| `EAS.1` | 236,513 | 0.04 |
| `EAS.2` | 1,272 | 0.01 |
| `n1` | 254,883 | 0.04 |
| `n2` | 17,694 | 0.00 |
| `root` | 22,906 | 0.00 |

log-Ne random-walk step scale tau_snp = 0.819

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,004.1 | 222 | 332.90 |
| SNP | +41.4 | 6 | 0.73 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.00 | -0.10 | +0.11 |
| **IBS** | -0.10 | +0.09 | +0.12 |
| **TSI** | +0.11 | +0.12 | -0.33 |

![spectrum](spectrum_fit.png)
