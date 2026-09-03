# `(((EAS,TSI.1),IBS),TSI.2)`

**Poisson, shared Ne** | topology 17 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5114.21 | +- 0.62 (MC) |
| logZ (importance sampling) | -5084.58 | |
| ESS of the IS weights | 10.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 13 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 1.1 +- 0.0 | 1.1 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 288.8 +- 0.4 | 289.9 |
| 3 | MERGE | n1 + IBS -> n2 | 12.6 +- 0.0 | 302.5 |
| 4 | MERGE | TSI.1 + n2 -> root | 1.0 +- 0.0 | 303.5 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,727,827 | 0.02 |
| `IBS` | 236,929 | 0.02 |
| `TSI` | 345,274 | 0.12 |
| `TSI.1` | 323,407 | 0.12 |
| `TSI.2` | 99 | 0.01 |
| `n1` | 70 | 0.01 |
| `n2` | 1,169 | 0.00 |
| `root` | 789 | 0.00 |

log-Ne random-walk step scale tau = 1.612

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,914.6 | 222 | 22994.55 |
| SNP | +38.9 | 6 | 1.58 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.30 | -0.49 | -0.10 |
| **IBS** | -0.49 | +0.14 | +0.85 |
| **TSI** | -0.10 | +0.85 | -0.60 |

![spectrum](spectrum_fit.png)
