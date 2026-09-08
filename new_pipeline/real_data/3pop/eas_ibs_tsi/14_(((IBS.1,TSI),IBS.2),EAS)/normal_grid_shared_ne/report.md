# `(((IBS.1,TSI),IBS.2),EAS)`

**Normal, shared Ne, recent grid** | topology 14 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +383.91 | +- 0.54 (MC) |
| logZ (importance sampling) | +417.51 | |
| ESS of the IS weights | 6.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 1 / 11 | 27 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 62.3 +- 0.2 | 62.3 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 71.8 +- 0.5 | 134.0 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 211.8 +- 0.9 | 345.9 |
| 4 | MERGE | n2 + EAS -> root | 11.5 +- 0.0 | 357.4 |

## Admixture fraction

**f = 0.117 +- 0.001** (fraction from `IBS.1`; 0.883 from `IBS.2`)

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 42,067,106 | 25,557,099 |
| `IBS` | 1,805,595 | 1,468,358 |
| `TSI` | 1,684,709 | 984,039 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 869,001 | 0.03 |
| `IBS` | 957,753 | 0.07 |
| `TSI` | 394,341 | 0.05 |
| `IBS.1` | 1,530 | 0.05 |
| `IBS.2` | 629,644 | 0.06 |
| `n1` | 44,515 | 0.05 |
| `n2` | 66 | 0.01 |
| `root` | 34 | 0.01 |

log-Ne random-walk step scale tau = 2.267

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +567.1 | 222 | 34.25 |
| SNP | +38.3 | 6 | 1.76 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.05 | +0.65 | -0.76 |
| **IBS** | +0.65 | -1.80 | +0.61 |
| **TSI** | -0.76 | +0.61 | +0.87 |

![spectrum](spectrum_fit.png)
