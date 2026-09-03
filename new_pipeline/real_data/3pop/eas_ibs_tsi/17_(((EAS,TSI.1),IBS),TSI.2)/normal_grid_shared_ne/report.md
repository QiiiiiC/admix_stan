# `(((EAS,TSI.1),IBS),TSI.2)`

**Normal, shared Ne, recent grid** | topology 17 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -36560.18 | +- 0.41 (MC) |
| logZ (importance sampling) | -36537.15 | |
| ESS of the IS weights | 5.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 1 | 24 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 196.2 +- 0.6 | 196.2 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 1.0 +- 0.0 | 197.2 |
| 3 | MERGE | n1 + IBS -> n2 | 1.8 +- 0.0 | 199.0 |
| 4 | MERGE | TSI.1 + n2 -> root | 1.1 +- 0.0 | 200.1 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,271,325 | 4,646,317 |
| `IBS` | 12,152,554 | 2,957,683 |
| `TSI` | 1,630,293 | 1,326,401 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,416,888 | 0.02 |
| `IBS` | 212,899 | 0.06 |
| `TSI` | 534,911 | 0.08 |
| `TSI.1` | 162 | 0.03 |
| `TSI.2` | 298 | 0.06 |
| `n1` | 298 | 0.06 |
| `n2` | 32,327 | 0.02 |
| `root` | 10,528 | 0.01 |

log-Ne random-walk step scale tau = 2.232

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -2,523.9 | 222 | 61.90 |
| SNP | -33,651.3 | 6 | 11231.62 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +109.50 | -119.02 | -97.33 |
| **IBS** | -119.02 | +93.12 | +143.14 |
| **TSI** | -97.33 | +143.14 | +50.93 |

![spectrum](spectrum_fit.png)
