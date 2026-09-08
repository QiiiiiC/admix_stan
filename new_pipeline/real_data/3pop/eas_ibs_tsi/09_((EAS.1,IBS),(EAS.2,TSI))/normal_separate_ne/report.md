# `((EAS.1,IBS),(EAS.2,TSI))`

**Normal, separate IBD/SNP Ne** | topology 09 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +515.31 | +- 0.29 (MC) |
| logZ (importance sampling) | +536.12 | |
| ESS of the IS weights | 3.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 2 / 11 | 21 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 120.6 +- 0.5 | 120.6 |
| 2 | MERGE | EAS.1 + IBS -> n1 | 4.3 +- 0.1 | 124.9 |
| 3 | MERGE | EAS.2 + TSI -> n2 | 18.8 +- 0.3 | 143.7 |
| 4 | MERGE | n1 + n2 -> root | 1.0 +- 0.0 | 144.7 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `EAS.1`; 1.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,751,950 | 0.07 |
| `IBS` | 654,752 | 0.05 |
| `TSI` | 330,149 | 0.03 |
| `EAS.1` | 7,019 | 0.03 |
| `EAS.2` | 41,830 | 0.02 |
| `n1` | 7,556 | 0.03 |
| `n2` | 224,561 | 0.03 |
| `root` | 83,285 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.949

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 673 | 0.01 |
| `IBS` | 1,498,779 | 0.12 |
| `TSI` | 73,520 | 0.15 |
| `EAS.1` | 34,534 | 0.01 |
| `EAS.2` | 20,064 | 0.01 |
| `n1` | 33,624 | 0.01 |
| `n2` | 27,466 | 0.01 |
| `root` | 28,247 | 0.01 |

log-Ne random-walk step scale tau_snp = 1.019

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +687.3 | 222 | 32.13 |
| SNP | +38.4 | 6 | 1.76 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.00 | -0.55 | +0.56 |
| **IBS** | -0.55 | +1.06 | -0.01 |
| **TSI** | +0.56 | -0.01 | -1.06 |

![spectrum](spectrum_fit.png)
