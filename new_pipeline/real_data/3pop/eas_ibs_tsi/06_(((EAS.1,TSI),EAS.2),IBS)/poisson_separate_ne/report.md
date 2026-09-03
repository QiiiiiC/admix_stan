# `(((EAS.1,TSI),EAS.2),IBS)`

**Poisson, separate IBD/SNP Ne** | topology 06 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -7842.73 | +- 0.89 (MC) |
| logZ (importance sampling) | -7796.01 | |
| ESS of the IS weights | 1.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 7 | 16 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 233.7 +- 0.3 | 234.7 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 1.5 +- 0.1 | 236.2 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 237.2 |

## Admixture fraction

**f = 0.999 +- 0.000** (fraction from `EAS.1`; 0.001 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,312,881 | 0.02 |
| `IBS` | 243,300 | 0.02 |
| `TSI` | 398,531 | 0.12 |
| `EAS.1` | 1,165,725 | 0.02 |
| `EAS.2` | 904 | 0.11 |
| `n1` | 104 | 0.11 |
| `n2` | 366 | 0.03 |
| `root` | 4,698 | 0.03 |

log-Ne random-walk step scale tau_ibd = 2.878

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,312 | 0.02 |
| `IBS` | 256,886 | 0.04 |
| `TSI` | 142,122 | 0.05 |
| `EAS.1` | 1,309 | 0.02 |
| `EAS.2` | 25,930 | 0.01 |
| `n1` | 26,001 | 0.01 |
| `n2` | 25,753 | 0.01 |
| `root` | 28,116 | 0.01 |

log-Ne random-walk step scale tau_snp = 1.012

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,583.3 | 222 | 307.66 |
| SNP | +27.5 | 6 | 5.38 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.00 | -0.29 | +0.30 |
| **IBS** | -0.29 | +0.61 | -0.06 |
| **TSI** | +0.30 | -0.06 | -0.53 |

![spectrum](spectrum_fit.png)
