# `((EAS,IBS.1),(IBS.2,TSI))`

**Normal, shared Ne, recent grid** | topology 15 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -42006.21 | +- 3.89 (MC) |
| logZ (importance sampling) | -41824.30 | |
| ESS of the IS weights | 1.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 7 | 16 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | IBS.1 + EAS -> n1 | 1.0 +- 0.0 | 12.0 |
| 3 | MERGE | IBS.2 + TSI -> n2 | 1.0 +- 0.0 | 13.0 |
| 4 | MERGE | n1 + n2 -> root | 364.5 +- 0.8 | 377.5 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `IBS.1`; 1.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 165,528,097 | 108,000,055 |
| `IBS` | 23,460,961 | 12,400,835 |
| `TSI` | 16,473,091 | 10,338,062 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 13,147,684 | 0.11 |
| `IBS` | 1,355,666 | 0.16 |
| `TSI` | 1,585,391 | 0.16 |
| `IBS.1` | 845,223 | 0.09 |
| `IBS.2` | 835,788 | 0.15 |
| `n1` | 845,122 | 0.09 |
| `n2` | 324,096 | 0.15 |
| `root` | 2 | 0.07 |

log-Ne random-walk step scale tau = 2.379

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -2,335.4 | 222 | 57.63 |
| SNP | -39,350.1 | 6 | 13131.23 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +116.63 | -114.79 | -115.73 |
| **IBS** | -114.79 | +112.96 | +113.25 |
| **TSI** | -115.73 | +113.25 | +114.14 |

![spectrum](spectrum_fit.png)
