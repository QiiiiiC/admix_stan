# `(((EAS,TSI),IBS.1),IBS.2)`

**Poisson, shared Ne, recent grid** | topology 12 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -43042.72 | +- 0.48 (MC) |
| logZ (importance sampling) | -43015.62 | |
| ESS of the IS weights | 9.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 1 | 9 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 214.9 +- 1.3 | 214.9 |
| 2 | MERGE | EAS + TSI -> n1 | 1.0 +- 0.0 | 215.9 |
| 3 | MERGE | IBS.2 + n1 -> n2 | 1.0 +- 0.0 | 216.9 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.0 +- 0.0 | 217.9 |

## Admixture fraction

**f = 0.001 +- 0.000** (fraction from `IBS.1`; 0.999 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 86,701,921 | 41,971,163 |
| `IBS` | 6,212,516 | 3,732,887 |
| `TSI` | 4,065,647 | 2,485,267 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,042,885 | 0.02 |
| `IBS` | 393,772 | 0.05 |
| `TSI` | 311,971 | 0.06 |
| `IBS.1` | 13,526 | 0.03 |
| `IBS.2` | 71 | 0.05 |
| `n1` | 1,025 | 0.04 |
| `n2` | 1,851 | 0.04 |
| `root` | 12,585 | 0.03 |

log-Ne random-walk step scale tau = 3.498

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -6,678.2 | 222 | 122.31 |
| SNP | -36,100.7 | 6 | 12048.10 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +112.61 | -96.94 | -125.75 |
| **IBS** | -96.94 | +41.20 | +153.71 |
| **TSI** | -125.75 | +153.71 | +95.16 |

![spectrum](spectrum_fit.png)
