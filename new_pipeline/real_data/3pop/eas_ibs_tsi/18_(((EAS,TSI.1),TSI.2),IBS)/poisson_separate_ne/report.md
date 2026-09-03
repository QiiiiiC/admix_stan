# `(((EAS,TSI.1),TSI.2),IBS)`

**Poisson, separate IBD/SNP Ne** | topology 18 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -4298.06 | +- 0.77 (MC) |
| logZ (importance sampling) | -4254.89 | |
| ESS of the IS weights | 2.3 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 13 | 16 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 131.5 +- 1.2 | 131.5 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 1.2 +- 0.0 | 132.6 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 44.1 +- 2.0 | 176.7 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 177.7 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,586,139 | 0.05 |
| `IBS` | 235,539 | 0.04 |
| `TSI` | 410,969 | 0.05 |
| `TSI.1` | 34,214 | 0.06 |
| `TSI.2` | 32,952 | 0.06 |
| `n1` | 32,952 | 0.06 |
| `n2` | 86,604 | 0.05 |
| `root` | 75,735 | 0.05 |

log-Ne random-walk step scale tau_ibd = 0.875

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 767 | 0.02 |
| `IBS` | 175,215 | 0.05 |
| `TSI` | 254,756 | 0.08 |
| `TSI.1` | 34,806 | 0.06 |
| `TSI.2` | 5,924 | 0.03 |
| `n1` | 5,923 | 0.03 |
| `n2` | 17,630 | 0.01 |
| `root` | 17,979 | 0.01 |

log-Ne random-walk step scale tau_snp = 1.251

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,120.2 | 222 | 77.79 |
| SNP | +31.9 | 6 | 3.90 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.11 | -0.19 | +0.40 |
| **IBS** | -0.19 | +0.21 | +0.16 |
| **TSI** | +0.40 | +0.16 | -0.92 |

![spectrum](spectrum_fit.png)
