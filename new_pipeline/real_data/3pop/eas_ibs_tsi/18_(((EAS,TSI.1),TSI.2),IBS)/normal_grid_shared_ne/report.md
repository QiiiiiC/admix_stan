# `(((EAS,TSI.1),TSI.2),IBS)`

**Normal, shared Ne, recent grid** | topology 18 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -35538.05 | +- 0.52 (MC) |
| logZ (importance sampling) | -35496.60 | |
| ESS of the IS weights | 1.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 1 / 9 | 27 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 43.6 +- 0.2 | 43.6 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 125.2 +- 0.2 | 168.8 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 2.8 +- 0.0 | 171.6 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 172.6 |

## Admixture fraction

**f = 0.998 +- 0.000** (fraction from `TSI.1`; 0.002 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 5,832,076 | 5,242,195 |
| `IBS` | 3,416,715 | 2,188,246 |
| `TSI` | 2,705,788 | 915,710 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,363,318 | 0.05 |
| `IBS` | 210,640 | 0.04 |
| `TSI` | 591,089 | 0.05 |
| `TSI.1` | 220,456 | 0.03 |
| `TSI.2` | 0 | 0.05 |
| `n1` | 1,404 | 0.02 |
| `n2` | 51,837 | 0.01 |
| `root` | 31,924 | 0.00 |

log-Ne random-walk step scale tau = 2.102

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -55.1 | 222 | 39.34 |
| SNP | -35,225.0 | 6 | 11756.21 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +111.42 | -123.36 | -96.78 |
| **IBS** | -123.36 | +92.86 | +152.29 |
| **TSI** | -96.78 | +152.29 | +41.26 |

![spectrum](spectrum_fit.png)
