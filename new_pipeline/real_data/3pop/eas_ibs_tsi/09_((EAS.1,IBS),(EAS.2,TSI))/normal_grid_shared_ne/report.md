# `((EAS.1,IBS),(EAS.2,TSI))`

**Normal, shared Ne, recent grid** | topology 09 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -36562.20 | +- 0.34 (MC) |
| logZ (importance sampling) | -36540.63 | |
| ESS of the IS weights | 11.8 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 1 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 194.9 +- 0.4 | 194.9 |
| 2 | MERGE | EAS.1 + IBS -> n1 | 1.6 +- 0.0 | 196.5 |
| 3 | MERGE | EAS.2 + TSI -> n2 | 1.0 +- 0.0 | 197.5 |
| 4 | MERGE | n1 + n2 -> root | 1.0 +- 0.0 | 198.5 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 12,522,802 | 3,512,709 |
| `IBS` | 10,879,693 | 6,209,034 |
| `TSI` | 248,321 | 166,849 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,405,146 | 0.03 |
| `IBS` | 212,496 | 0.02 |
| `TSI` | 668,372 | 0.07 |
| `EAS.1` | 290 | 0.04 |
| `EAS.2` | 42 | 0.06 |
| `n1` | 90,522 | 0.02 |
| `n2` | 42 | 0.06 |
| `root` | 11,177 | 0.02 |

log-Ne random-walk step scale tau = 4.839

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -2,241.1 | 222 | 59.22 |
| SNP | -33,911.0 | 6 | 11318.22 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +109.91 | -119.49 | -97.67 |
| **IBS** | -119.49 | +93.48 | +143.70 |
| **TSI** | -97.67 | +143.70 | +51.05 |

![spectrum](spectrum_fit.png)
