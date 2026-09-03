# `(((EAS.1,TSI),EAS.2),IBS)`

**Normal, shared Ne** | topology 06 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -39408.11 | +- 1.99 (MC) |
| logZ (importance sampling) | -39289.33 | |
| ESS of the IS weights | 1.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 7 | 18 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 163.2 +- 6.5 | 163.2 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 6.1 +- 3.2 | 169.2 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 1.2 +- 0.2 | 170.5 |
| 4 | MERGE | n2 + IBS -> root | 1.1 +- 0.0 | 171.5 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `EAS.1`; 1.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,526,416 | 0.10 |
| `IBS` | 210,145 | 0.05 |
| `TSI` | 320,730 | 0.11 |
| `EAS.1` | 107,621 | 0.44 |
| `EAS.2` | 3,865 | 0.76 |
| `n1` | 28,636 | 0.26 |
| `n2` | 30,997 | 0.24 |
| `root` | 33,745 | 0.11 |

log-Ne random-walk step scale tau = 0.974

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -315.5 | 222 | 41.42 |
| SNP | -38,880.8 | 6 | 12974.79 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +116.23 | -114.23 | -115.51 |
| **IBS** | -114.23 | +110.30 | +114.94 |
| **TSI** | -115.51 | +114.94 | +112.12 |

![spectrum](spectrum_fit.png)
