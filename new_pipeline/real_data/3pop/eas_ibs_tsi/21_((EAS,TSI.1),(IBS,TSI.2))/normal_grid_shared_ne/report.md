# `((EAS,TSI.1),(IBS,TSI.2))`

**Normal, shared Ne, recent grid** | topology 21 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -42059.23 | +- 1.31 (MC) |
| logZ (importance sampling) | -41977.76 | |
| ESS of the IS weights | 1.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 13 | 18 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | TSI.1 + EAS -> n1 | 1.0 +- 0.0 | 12.0 |
| 3 | MERGE | TSI.2 + IBS -> n2 | 1.0 +- 0.0 | 13.0 |
| 4 | MERGE | n1 + n2 -> root | 361.4 +- 1.7 | 374.4 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `TSI.1`; 1.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 49,338,297 | 29,320,994 |
| `IBS` | 903,029 | 779,015 |
| `TSI` | 1,048,952 | 911,292 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 3,389,092 | 0.11 |
| `IBS` | 464,984 | 0.08 |
| `TSI` | 515,320 | 0.08 |
| `TSI.1` | 861,879 | 0.05 |
| `TSI.2` | 452,086 | 0.08 |
| `n1` | 860,783 | 0.05 |
| `n2` | 328,930 | 0.09 |
| `root` | 4 | 0.11 |

log-Ne random-walk step scale tau = 1.788

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -2,267.4 | 222 | 56.86 |
| SNP | -39,367.8 | 6 | 13137.13 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +116.66 | -114.83 | -115.77 |
| **IBS** | -114.83 | +112.95 | +113.32 |
| **TSI** | -115.77 | +113.32 | +114.13 |

![spectrum](spectrum_fit.png)
