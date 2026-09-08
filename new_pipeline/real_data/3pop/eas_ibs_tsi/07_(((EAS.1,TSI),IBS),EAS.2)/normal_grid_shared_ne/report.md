# `(((EAS.1,TSI),IBS),EAS.2)`

**Normal, shared Ne, recent grid** | topology 07 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +3813.59 | +- 0.12 (MC) |
| logZ (importance sampling) | +3824.37 | |
| ESS of the IS weights | 4.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 1 / 7 | 107 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 134.4 +- 0.4 | 134.4 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 20.8 +- 0.5 | 155.2 |
| 3 | MERGE | n1 + IBS -> n2 | 1.4 +- 0.0 | 156.6 |
| 4 | MERGE | EAS.1 + n2 -> root | 1,975.6 +- 10.0 | 2,132.2 |

## Admixture fraction

**f = 0.987 +- 0.000** (fraction from `EAS.1`; 0.013 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,997,235 | 4,474,199 |
| `IBS` | 2,568,215 | 1,578,565 |
| `TSI` | 675,129 | 561,407 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,348,671 | 0.01 |
| `IBS` | 275,953 | 0.02 |
| `TSI` | 416,548 | 0.04 |
| `EAS.1` | 35,706 | 0.02 |
| `EAS.2` | 18,579 | 0.01 |
| `n1` | 19,396 | 0.01 |
| `n2` | 15,179 | 0.01 |
| `root` | 14,037 | 0.00 |

log-Ne random-walk step scale tau = 1.560

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +3,958.5 | 222 | 3.21 |
| SNP | +33.1 | 6 | 3.51 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.34 | +0.44 | +0.22 |
| **IBS** | +0.44 | +1.86 | -2.89 |
| **TSI** | +0.22 | -2.89 | +2.30 |

![spectrum](spectrum_fit.png)
