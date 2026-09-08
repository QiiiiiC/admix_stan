# `(((EAS,TSI),IBS.1),IBS.2)`

**Poisson, separate IBD/SNP Ne** | topology 12 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -7156.30 | +- 0.46 (MC) |
| logZ (importance sampling) | -7122.52 | |
| ESS of the IS weights | 2.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 1 / 3 | 19 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 82.5 +- 0.9 | 82.5 |
| 2 | MERGE | EAS + TSI -> n1 | 144.9 +- 0.6 | 227.4 |
| 3 | MERGE | IBS.2 + n1 -> n2 | 1.2 +- 0.0 | 228.6 |
| 4 | MERGE | IBS.1 + n2 -> root | 342.2 +- 10.4 | 570.8 |

## Admixture fraction

**f = 0.419 +- 0.006** (fraction from `IBS.1`; 0.581 from `IBS.2`)

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,203,758 | 0.02 |
| `IBS` | 857,989 | 0.04 |
| `TSI` | 344,150 | 0.03 |
| `IBS.1` | 52,244 | 0.04 |
| `IBS.2` | 19,415 | 0.03 |
| `n1` | 4,333 | 0.05 |
| `n2` | 4,148 | 0.05 |
| `root` | 9,612 | 0.02 |

log-Ne random-walk step scale tau_ibd = 1.270

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,261 | 0.02 |
| `IBS` | 828,792 | 0.03 |
| `TSI` | 177,069 | 0.05 |
| `IBS.1` | 265,140 | 0.04 |
| `IBS.2` | 344,829 | 0.04 |
| `n1` | 54,650 | 0.03 |
| `n2` | 59,800 | 0.03 |
| `root` | 20,655 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.973

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -6,983.8 | 222 | 225.75 |
| SNP | +30.7 | 6 | 4.32 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.02 | +0.23 | -0.19 |
| **IBS** | +0.23 | -0.65 | +0.23 |
| **TSI** | -0.19 | +0.23 | +0.15 |

![spectrum](spectrum_fit.png)
