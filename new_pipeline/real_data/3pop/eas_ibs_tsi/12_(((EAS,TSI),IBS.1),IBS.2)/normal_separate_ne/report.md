# `(((EAS,TSI),IBS.1),IBS.2)`

**Normal, separate IBD/SNP Ne** | topology 12 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1131.40 | +- 0.42 (MC) |
| logZ (importance sampling) | -1102.51 | |
| ESS of the IS weights | 6.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 1 / 11 | 21 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 84.5 +- 2.4 | 84.5 |
| 2 | MERGE | EAS + TSI -> n1 | 74.1 +- 4.4 | 158.7 |
| 3 | MERGE | IBS.2 + n1 -> n2 | 1.5 +- 0.3 | 160.1 |
| 4 | MERGE | IBS.1 + n2 -> root | 73.3 +- 59.5 | 233.5 |

## Admixture fraction

**f = 0.306 +- 0.006** (fraction from `IBS.1`; 0.694 from `IBS.2`)

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,575,090 | 0.05 |
| `IBS` | 834,373 | 0.13 |
| `TSI` | 367,311 | 0.03 |
| `IBS.1` | 52,366 | 0.55 |
| `IBS.2` | 24,708 | 0.15 |
| `n1` | 23,251 | 0.27 |
| `n2` | 21,494 | 0.23 |
| `root` | 66,675 | 0.05 |

log-Ne random-walk step scale tau_ibd = 1.282

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 880 | 0.02 |
| `IBS` | 167,161 | 0.04 |
| `TSI` | 123,564 | 0.12 |
| `IBS.1` | 47,962 | 0.10 |
| `IBS.2` | 73,394 | 0.11 |
| `n1` | 31,599 | 0.10 |
| `n2` | 31,329 | 0.09 |
| `root` | 25,088 | 0.25 |

log-Ne random-walk step scale tau_snp = 1.042

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -979.1 | 222 | 47.26 |
| SNP | +32.8 | 6 | 3.60 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.05 | +0.14 | -0.24 |
| **IBS** | +0.14 | -0.54 | +0.29 |
| **TSI** | -0.24 | +0.29 | +0.19 |

![spectrum](spectrum_fit.png)
