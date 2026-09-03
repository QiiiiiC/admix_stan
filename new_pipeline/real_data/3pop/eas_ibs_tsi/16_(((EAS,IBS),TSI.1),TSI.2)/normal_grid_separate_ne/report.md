# `(((EAS,IBS),TSI.1),TSI.2)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 16 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1406.52 | +- 0.59 (MC) |
| logZ (importance sampling) | -1363.09 | |
| ESS of the IS weights | 3.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 13 | 25 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 165.4 +- 1.5 | 165.4 |
| 2 | MERGE | EAS + IBS -> n1 | 3.3 +- 0.7 | 168.7 |
| 3 | MERGE | TSI.2 + n1 -> n2 | 1.8 +- 0.0 | 170.5 |
| 4 | MERGE | TSI.1 + n2 -> root | 5.7 +- 0.6 | 176.2 |

## Admixture fraction

**f = 0.087 +- 0.002** (fraction from `TSI.1`; 0.913 from `TSI.2`)

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 24,295,646 | 17,583,233 |
| `IBS` | 4,311,913 | 2,663,950 |
| `TSI` | 2,433,010 | 1,439,408 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,414,663 | 0.03 |
| `IBS` | 234,079 | 0.06 |
| `TSI` | 376,272 | 0.06 |
| `TSI.1` | 21,302 | 0.16 |
| `TSI.2` | 3,849 | 0.09 |
| `n1` | 5,328 | 0.07 |
| `n2` | 5,527 | 0.04 |
| `root` | 85,714 | 0.26 |

log-Ne random-walk step scale tau_ibd = 2.303

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,144 | 1,129 |
| `IBS` | 119,472 | 118,042 |
| `TSI` | 200,310 | 183,885 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 925 | 0.02 |
| `IBS` | 106,993 | 0.06 |
| `TSI` | 174,189 | 0.17 |
| `TSI.1` | 14,361 | 0.06 |
| `TSI.2` | 19,849 | 0.05 |
| `n1` | 18,898 | 0.06 |
| `n2` | 17,216 | 0.06 |
| `root` | 20,167 | 0.07 |

log-Ne random-walk step scale tau_snp = 0.588

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,127.3 | 222 | 48.69 |
| SNP | +35.4 | 6 | 2.73 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.00 | +0.19 | -0.20 |
| **IBS** | +0.19 | -0.58 | +0.22 |
| **TSI** | -0.20 | +0.22 | +0.18 |

![spectrum](spectrum_fit.png)
