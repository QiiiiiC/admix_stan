# `(((EAS.1,TSI),IBS),EAS.2)`

**Poisson, shared Ne** | topology 07 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -6052.45 | +- 1.65 (MC) |
| logZ (importance sampling) | -5963.35 | |
| ESS of the IS weights | 3.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 7 | 18 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 255.5 +- 1.7 | 256.5 |
| 3 | MERGE | n1 + IBS -> n2 | 1.0 +- 0.0 | 257.5 |
| 4 | MERGE | EAS.1 + n2 -> root | 67.6 +- 3.0 | 325.2 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 861,196 | 0.06 |
| `IBS` | 317,875 | 0.09 |
| `TSI` | 420,753 | 0.13 |
| `EAS.1` | 888,765 | 0.06 |
| `EAS.2` | 350 | 0.25 |
| `n1` | 489 | 0.05 |
| `n2` | 377 | 0.05 |
| `root` | 591 | 0.12 |

log-Ne random-walk step scale tau = 1.168

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,869.7 | 222 | 185.24 |
| SNP | +4.3 | 6 | 13.10 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.12 | -0.95 | +0.72 |
| **IBS** | -0.95 | +0.62 | +1.28 |
| **TSI** | +0.72 | +1.28 | -2.57 |

![spectrum](spectrum_fit.png)
