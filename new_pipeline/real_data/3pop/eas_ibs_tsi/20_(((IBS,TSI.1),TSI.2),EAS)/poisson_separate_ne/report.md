# `(((IBS,TSI.1),TSI.2),EAS)`

**Poisson, separate IBD/SNP Ne** | topology 20 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5526.59 | +- 0.42 (MC) |
| logZ (importance sampling) | -5496.62 | |
| ESS of the IS weights | 2.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 2 / 7 | 18 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 16.0 +- 0.7 | 16.0 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 51.7 +- 2.2 | 67.7 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 80.2 +- 1.0 | 148.0 |
| 4 | MERGE | n2 + EAS -> root | 146.5 +- 1.0 | 294.5 |

## Admixture fraction

**f = 0.999 +- 0.000** (fraction from `TSI.1`; 0.001 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 894,869 | 0.02 |
| `IBS` | 988,413 | 0.09 |
| `TSI` | 624,566 | 0.07 |
| `TSI.1` | 382,901 | 0.07 |
| `TSI.2` | 104,663 | 0.05 |
| `n1` | 82,456 | 0.05 |
| `n2` | 25,454 | 0.04 |
| `root` | 1,693 | 0.06 |

log-Ne random-walk step scale tau_ibd = 1.424

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,650 | 0.01 |
| `IBS` | 122,036 | 0.03 |
| `TSI` | 105,944 | 0.03 |
| `TSI.1` | 102,209 | 0.03 |
| `TSI.2` | 90,210 | 0.03 |
| `n1` | 97,297 | 0.03 |
| `n2` | 78,238 | 0.03 |
| `root` | 23,000 | 0.02 |

log-Ne random-walk step scale tau_snp = 0.497

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,360.0 | 222 | 43.77 |
| SNP | +38.2 | 6 | 1.80 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.02 | +0.06 | -0.01 |
| **IBS** | +0.06 | -0.23 | +0.13 |
| **TSI** | -0.01 | +0.13 | -0.10 |

![spectrum](spectrum_fit.png)
