# `(((EAS,IBS),TSI.1),TSI.2)`

**Poisson, shared Ne** | topology 16 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -43686.04 | +- 1.09 (MC) |
| logZ (importance sampling) | -43657.12 | |
| ESS of the IS weights | 6.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 7 | 18 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | EAS + IBS -> n1 | 282.8 +- 0.6 | 283.8 |
| 3 | MERGE | TSI.2 + n1 -> n2 | 7.1 +- 0.1 | 290.9 |
| 4 | MERGE | TSI.1 + n2 -> root | 1.0 +- 0.0 | 291.9 |

## Admixture fraction

**f = 0.135 +- 0.005** (fraction from `TSI.1`; 0.865 from `TSI.2`)

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,200,841 | 0.03 |
| `IBS` | 244,223 | 0.04 |
| `TSI` | 223,962 | 0.08 |
| `TSI.1` | 36,328 | 0.48 |
| `TSI.2` | 302,587 | 0.13 |
| `n1` | 251 | 0.02 |
| `n2` | 620 | 0.01 |
| `root` | 724 | 0.01 |

log-Ne random-walk step scale tau = 1.446

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,593.8 | 222 | 6361.72 |
| SNP | -35,900.3 | 6 | 11981.31 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +112.54 | -124.33 | -98.02 |
| **IBS** | -124.33 | +93.04 | +154.06 |
| **TSI** | -98.02 | +154.06 | +41.98 |

![spectrum](spectrum_fit.png)
