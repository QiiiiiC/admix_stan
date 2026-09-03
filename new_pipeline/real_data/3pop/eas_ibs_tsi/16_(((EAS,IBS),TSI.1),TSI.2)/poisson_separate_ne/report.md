# `(((EAS,IBS),TSI.1),TSI.2)`

**Poisson, separate IBD/SNP Ne** | topology 16 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -7694.85 | +- 0.39 (MC) |
| logZ (importance sampling) | -7668.79 | |
| ESS of the IS weights | 4.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 7 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | EAS + IBS -> n1 | 227.8 +- 1.1 | 228.8 |
| 3 | MERGE | TSI.2 + n1 -> n2 | 8.7 +- 0.7 | 237.6 |
| 4 | MERGE | TSI.1 + n2 -> root | 1.0 +- 0.0 | 238.6 |

## Admixture fraction

**f = 0.176 +- 0.001** (fraction from `TSI.1`; 0.824 from `TSI.2`)

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,228,666 | 0.01 |
| `IBS` | 254,344 | 0.05 |
| `TSI` | 212,584 | 0.04 |
| `TSI.1` | 37,188 | 0.02 |
| `TSI.2` | 308,791 | 0.05 |
| `n1` | 2,128 | 0.05 |
| `n2` | 4,596 | 0.01 |
| `root` | 5,161 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.440

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,269 | 0.02 |
| `IBS` | 171,340 | 0.17 |
| `TSI` | 172,358 | 0.08 |
| `TSI.1` | 54,004 | 0.03 |
| `TSI.2` | 221,826 | 0.10 |
| `n1` | 14,349 | 0.01 |
| `n2` | 18,035 | 0.01 |
| `root` | 18,579 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.805

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,530.0 | 222 | 345.59 |
| SNP | +34.0 | 6 | 3.22 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.03 | +0.04 | +0.03 |
| **IBS** | +0.04 | -0.25 | +0.18 |
| **TSI** | +0.03 | +0.18 | -0.22 |

![spectrum](spectrum_fit.png)
