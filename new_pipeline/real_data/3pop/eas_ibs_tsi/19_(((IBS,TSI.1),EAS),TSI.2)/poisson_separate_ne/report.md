# `(((IBS,TSI.1),EAS),TSI.2)`

**Poisson, separate IBD/SNP Ne** | topology 19 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5870.47 | +- 0.48 (MC) |
| logZ (importance sampling) | -5834.23 | |
| ESS of the IS weights | 7.3 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 13 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 159.8 +- 1.4 | 160.8 |
| 3 | MERGE | n1 + EAS -> n2 | 152.8 +- 9.9 | 313.6 |
| 4 | MERGE | TSI.1 + n2 -> root | 8.0 +- 0.7 | 321.6 |

## Admixture fraction

**f = 0.001 +- 0.000** (fraction from `TSI.1`; 0.999 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 895,714 | 0.02 |
| `IBS` | 313,306 | 0.07 |
| `TSI` | 426,506 | 0.04 |
| `TSI.1` | 36 | 0.24 |
| `TSI.2` | 424,110 | 0.04 |
| `n1` | 13,892 | 0.05 |
| `n2` | 1,131 | 0.43 |
| `root` | 822 | 0.46 |

log-Ne random-walk step scale tau_ibd = 1.283

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 28,504 | 0.17 |
| `IBS` | 103,986 | 0.09 |
| `TSI` | 105,072 | 0.16 |
| `TSI.1` | 1,046 | 0.34 |
| `TSI.2` | 104,381 | 0.16 |
| `n1` | 905 | 0.08 |
| `n2` | 9,025 | 0.03 |
| `root` | 8,648 | 0.05 |

log-Ne random-walk step scale tau_snp = 1.217

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,685.9 | 222 | 46.46 |
| SNP | +34.1 | 6 | 3.18 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.02 | +0.14 | -0.18 |
| **IBS** | +0.14 | -0.86 | +0.64 |
| **TSI** | -0.18 | +0.64 | -0.26 |

![spectrum](spectrum_fit.png)
