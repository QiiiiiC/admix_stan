# `(((IBS.1,TSI),EAS),IBS.2)`

**Normal, separate IBD/SNP Ne** | topology 13 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -83.94 | +- 0.88 (MC) |
| logZ (importance sampling) | -32.01 | |
| ESS of the IS weights | 3.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 7 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 87.7 +- 1.2 | 87.7 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 44.6 +- 0.4 | 132.2 |
| 3 | MERGE | n1 + EAS -> n2 | 236.9 +- 0.2 | 369.1 |
| 4 | MERGE | IBS.1 + n2 -> root | 208.6 +- 0.3 | 577.7 |

## Admixture fraction

**f = 0.015 +- 0.002** (fraction from `IBS.1`; 0.985 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 976,274 | 0.03 |
| `IBS` | 835,939 | 0.01 |
| `TSI` | 403,785 | 0.07 |
| `IBS.1` | 0 | 0.05 |
| `IBS.2` | 63,355 | 0.04 |
| `n1` | 52,685 | 0.04 |
| `n2` | 3 | 0.01 |
| `root` | 2,397 | 0.00 |

log-Ne random-walk step scale tau_ibd = 1.590

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 13,287,218 | 0.03 |
| `IBS` | 1,096,260 | 0.03 |
| `TSI` | 9,386,107 | 0.03 |
| `IBS.1` | 23,864 | 0.01 |
| `IBS.2` | 11,315 | 0.02 |
| `n1` | 1,301 | 0.03 |
| `n2` | 2,058 | 0.02 |
| `root` | 25,393 | 0.01 |

log-Ne random-walk step scale tau_snp = 1.697

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +168.6 | 222 | 37.65 |
| SNP | +11.8 | 6 | 10.59 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.34 | -1.00 | +0.34 |
| **IBS** | -1.00 | -0.16 | +2.21 |
| **TSI** | +0.34 | +2.21 | -2.72 |

![spectrum](spectrum_fit.png)
