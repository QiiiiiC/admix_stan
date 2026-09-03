# `(((IBS.1,TSI),EAS),IBS.2)`

**Poisson, separate IBD/SNP Ne** | topology 13 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5877.34 | +- 0.43 (MC) |
| logZ (importance sampling) | -5843.71 | |
| ESS of the IS weights | 2.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 7 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 159.6 +- 1.2 | 159.6 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 1.0 +- 0.0 | 160.6 |
| 3 | MERGE | n1 + EAS -> n2 | 147.2 +- 14.8 | 307.8 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.1 +- 0.0 | 308.9 |

## Admixture fraction

**f = 0.002 +- 0.000** (fraction from `IBS.1`; 0.998 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 895,904 | 0.02 |
| `IBS` | 318,650 | 0.05 |
| `TSI` | 424,966 | 0.07 |
| `IBS.1` | 1,077 | 0.61 |
| `IBS.2` | 13,291 | 0.07 |
| `n1` | 13,887 | 0.04 |
| `n2` | 1,289 | 0.63 |
| `root` | 1,172 | 0.62 |

log-Ne random-walk step scale tau_ibd = 1.373

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 493,613 | 0.83 |
| `IBS` | 338,466 | 0.57 |
| `TSI` | 121,538 | 0.08 |
| `IBS.1` | 9,798 | 0.26 |
| `IBS.2` | 884 | 0.10 |
| `n1` | 822 | 0.11 |
| `n2` | 10,742 | 0.23 |
| `root` | 10,321 | 0.25 |

log-Ne random-walk step scale tau_snp = 1.296

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,676.5 | 222 | 46.47 |
| SNP | +37.3 | 6 | 2.10 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.08 | -0.09 | -0.06 |
| **IBS** | -0.09 | -0.45 | +0.67 |
| **TSI** | -0.06 | +0.67 | -0.52 |

![spectrum](spectrum_fit.png)
