# `(((EAS.1,TSI),EAS.2),IBS)`

**Normal, shared Ne, recent grid** | topology 06 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -36501.95 | +- 0.21 (MC) |
| logZ (importance sampling) | -36483.31 | |
| ESS of the IS weights | 1.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 1 / 7 | 60 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 195.4 +- 0.3 | 195.4 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 1.1 +- 0.0 | 196.5 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 2.1 +- 0.0 | 198.5 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 199.5 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,905,176 | 3,080,930 |
| `IBS` | 7,944,026 | 4,168,320 |
| `TSI` | 124,960 | 143,833 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,444,456 | 0.03 |
| `IBS` | 212,524 | 0.03 |
| `TSI` | 729,223 | 0.06 |
| `EAS.1` | 539 | 0.03 |
| `EAS.2` | 70 | 0.03 |
| `n1` | 86 | 0.04 |
| `n2` | 24,271 | 0.03 |
| `root` | 10,683 | 0.01 |

log-Ne random-walk step scale tau = 4.066

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -2,421.8 | 222 | 60.81 |
| SNP | -33,773.4 | 6 | 11272.33 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +109.70 | -119.23 | -97.51 |
| **IBS** | -119.23 | +93.17 | +143.52 |
| **TSI** | -97.51 | +143.52 | +50.91 |

![spectrum](spectrum_fit.png)
