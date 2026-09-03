# `(((EAS,IBS.1),TSI),IBS.2)`

**Normal, separate IBD/SNP Ne** | topology 11 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -26035.18 | +- 0.27 (MC) |
| logZ (importance sampling) | -26015.45 | |
| ESS of the IS weights | 2.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 7 | 18 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 105.1 +- 0.7 | 105.1 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 28.4 +- 0.4 | 133.5 |
| 3 | MERGE | n1 + TSI -> n2 | 10.5 +- 0.1 | 144.0 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.0 +- 0.0 | 145.0 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `IBS.1`; 0.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,651,682 | 0.04 |
| `IBS` | 723,338 | 0.09 |
| `TSI` | 331,071 | 0.05 |
| `IBS.1` | 20,328 | 0.02 |
| `IBS.2` | 15,642 | 0.03 |
| `n1` | 15,643 | 0.03 |
| `n2` | 132,013 | 0.01 |
| `root` | 82,841 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.514

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 5,106 | 0.01 |
| `IBS` | 5,106 | 0.01 |
| `TSI` | 5,106 | 0.01 |
| `IBS.1` | 5,106 | 0.01 |
| `IBS.2` | 5,106 | 0.01 |
| `n1` | 5,106 | 0.01 |
| `n2` | 5,106 | 0.01 |
| `root` | 5,106 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.000

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +700.3 | 222 | 32.03 |
| SNP | -26,514.7 | 6 | 8852.77 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +90.11 | -88.48 | -89.62 |
| **IBS** | -88.48 | +9.85 | +170.06 |
| **TSI** | -89.62 | +170.06 | +10.90 |

![spectrum](spectrum_fit.png)
