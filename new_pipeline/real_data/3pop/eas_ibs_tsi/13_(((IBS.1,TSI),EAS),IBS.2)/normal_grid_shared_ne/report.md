# `(((IBS.1,TSI),EAS),IBS.2)`

**Normal, shared Ne, recent grid** | topology 13 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1099.10 | +- 0.63 (MC) |
| logZ (importance sampling) | -1058.32 | |
| ESS of the IS weights | 2.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 1 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 208.0 +- 0.7 | 208.0 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 3.1 +- 0.0 | 211.1 |
| 3 | MERGE | n1 + EAS -> n2 | 208.5 +- 1.0 | 419.7 |
| 4 | MERGE | IBS.1 + n2 -> root | 195.2 +- 0.3 | 614.9 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `IBS.1`; 1.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 127,784,131 | 49,397,667 |
| `IBS` | 2,558,018 | 1,733,137 |
| `TSI` | 551,345 | 474,108 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 717,764 | 0.03 |
| `IBS` | 375,001 | 0.02 |
| `TSI` | 478,032 | 0.03 |
| `IBS.1` | 2,899 | 0.01 |
| `IBS.2` | 1,539 | 0.02 |
| `n1` | 1,176 | 0.03 |
| `n2` | 1 | 0.04 |
| `root` | 7,661 | 0.01 |

log-Ne random-walk step scale tau = 2.206

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -872.4 | 222 | 47.64 |
| SNP | +12.8 | 6 | 10.28 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +1.59 | -0.52 | -2.63 |
| **IBS** | -0.52 | -0.95 | +2.08 |
| **TSI** | -2.63 | +2.08 | +3.05 |

![spectrum](spectrum_fit.png)
