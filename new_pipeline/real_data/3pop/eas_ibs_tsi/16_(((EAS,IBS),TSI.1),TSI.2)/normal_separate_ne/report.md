# `(((EAS,IBS),TSI.1),TSI.2)`

**Normal, separate IBD/SNP Ne** | topology 16 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1526.53 | +- 0.51 (MC) |
| logZ (importance sampling) | -1497.22 | |
| ESS of the IS weights | 2.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 2 / 6 | 25 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 133.4 +- 0.5 | 133.4 |
| 2 | MERGE | EAS + IBS -> n1 | 30.7 +- 0.1 | 164.1 |
| 3 | MERGE | TSI.2 + n1 -> n2 | 1.3 +- 0.0 | 165.4 |
| 4 | MERGE | TSI.1 + n2 -> root | 20.2 +- 0.2 | 185.6 |

## Admixture fraction

**f = 0.309 +- 0.007** (fraction from `TSI.1`; 0.691 from `TSI.2`)

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,616,788 | 0.04 |
| `IBS` | 245,931 | 0.03 |
| `TSI` | 419,610 | 0.09 |
| `TSI.1` | 246,684 | 0.05 |
| `TSI.2` | 15,042 | 0.05 |
| `n1` | 11,637 | 0.04 |
| `n2` | 11,817 | 0.04 |
| `root` | 130,819 | 0.04 |

log-Ne random-walk step scale tau_ibd = 1.341

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 910 | 0.01 |
| `IBS` | 121,107 | 0.05 |
| `TSI` | 162,524 | 0.03 |
| `TSI.1` | 21,493 | 0.04 |
| `TSI.2` | 39,807 | 0.01 |
| `n1` | 24,632 | 0.00 |
| `n2` | 24,536 | 0.00 |
| `root` | 24,824 | 0.00 |

log-Ne random-walk step scale tau_snp = 0.958

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,369.4 | 222 | 50.65 |
| SNP | +40.7 | 6 | 0.98 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.03 | -0.07 | +0.01 |
| **IBS** | -0.07 | -0.19 | +0.34 |
| **TSI** | +0.01 | +0.34 | -0.34 |

![spectrum](spectrum_fit.png)
