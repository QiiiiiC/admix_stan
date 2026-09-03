# `(((EAS,IBS.1),TSI),IBS.2)`

**Poisson, shared Ne, recent grid** | topology 11 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5019.70 | +- 1.16 (MC) |
| logZ (importance sampling) | -4951.79 | |
| ESS of the IS weights | 5.3 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 1 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 276.9 +- 3.6 | 287.9 |
| 3 | MERGE | n1 + TSI -> n2 | 19.1 +- 7.8 | 307.0 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.0 +- 0.0 | 308.0 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `IBS.1`; 0.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 5,474,401 | 4,386,372 |
| `IBS` | 730,937 | 578,034 |
| `TSI` | 620,899 | 535,159 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,627,613 | 0.03 |
| `IBS` | 206,828 | 0.10 |
| `TSI` | 318,492 | 0.09 |
| `IBS.1` | 229,722 | 0.10 |
| `IBS.2` | 106 | 0.50 |
| `n1` | 106 | 0.50 |
| `n2` | 739 | 0.12 |
| `root` | 670 | 0.16 |

log-Ne random-walk step scale tau = 1.201

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,757.3 | 222 | 30768.31 |
| SNP | +26.5 | 6 | 5.69 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.27 | +0.63 | -1.17 |
| **IBS** | +0.63 | -2.75 | +1.66 |
| **TSI** | -1.17 | +1.66 | +0.67 |

![spectrum](spectrum_fit.png)
