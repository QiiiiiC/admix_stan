# `(((EAS.1,IBS),TSI),EAS.2)`

**Normal, separate IBD/SNP Ne** | topology 05 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +696.10 | +- 0.36 (MC) |
| logZ (importance sampling) | +728.55 | |
| ESS of the IS weights | 1.3 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 13 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 63.2 +- 0.6 | 63.2 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 1.0 +- 0.0 | 64.2 |
| 3 | MERGE | n1 + TSI -> n2 | 1.0 +- 0.0 | 65.2 |
| 4 | MERGE | EAS.1 + n2 -> root | 375.3 +- 1.5 | 440.5 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 4,248,716 | 0.05 |
| `IBS` | 245,762 | 0.05 |
| `TSI` | 382,120 | 0.05 |
| `EAS.1` | 271,353 | 0.02 |
| `EAS.2` | 179,149 | 0.04 |
| `n1` | 179,201 | 0.04 |
| `n2` | 214,121 | 0.03 |
| `root` | 2 | 0.20 |

log-Ne random-walk step scale tau_ibd = 1.717

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 5,107 | 0.01 |
| `IBS` | 5,106 | 0.01 |
| `TSI` | 5,104 | 0.01 |
| `EAS.1` | 5,108 | 0.01 |
| `EAS.2` | 5,096 | 0.01 |
| `n1` | 5,096 | 0.01 |
| `n2` | 5,095 | 0.01 |
| `root` | 5,191 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.006

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +2,512.3 | 222 | 15.48 |
| SNP | -1,532.0 | 6 | 525.21 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +9.60 | -9.39 | -9.58 |
| **IBS** | -9.39 | -22.43 | +43.21 |
| **TSI** | -9.58 | +43.21 | -22.39 |

![spectrum](spectrum_fit.png)
