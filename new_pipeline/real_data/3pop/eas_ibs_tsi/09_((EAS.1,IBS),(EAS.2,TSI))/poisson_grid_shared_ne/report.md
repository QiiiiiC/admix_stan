# `((EAS.1,IBS),(EAS.2,TSI))`

**Poisson, shared Ne, recent grid** | topology 09 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -4409.02 | +- 3.18 (MC) |
| logZ (importance sampling) | -4238.51 | |
| ESS of the IS weights | 7.4 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 1 / 4 | 28 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 146.7 +- 2.3 | 146.7 |
| 2 | MERGE | EAS.1 + IBS -> n1 | 59.3 +- 5.8 | 205.9 |
| 3 | MERGE | EAS.2 + TSI -> n2 | 23.7 +- 4.0 | 229.6 |
| 4 | MERGE | n1 + n2 -> root | 4.9 +- 0.4 | 234.4 |

## Admixture fraction

**f = 0.967 +- 0.001** (fraction from `EAS.1`; 0.033 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 3,102,275 | 4,549,983 |
| `IBS` | 3,987,277 | 2,814,994 |
| `TSI` | 8,369,368 | 3,473,252 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,427,896 | 0.10 |
| `IBS` | 224,467 | 0.10 |
| `TSI` | 306,937 | 0.13 |
| `EAS.1` | 947,962 | 0.24 |
| `EAS.2` | 0 | 0.14 |
| `n1` | 2,480,285 | 0.43 |
| `n2` | 22,609 | 0.40 |
| `root` | 10,633 | 0.39 |

log-Ne random-walk step scale tau = 2.282

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,172.3 | 222 | 527.30 |
| SNP | -0.6 | 6 | 14.73 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.13 | +0.15 | +0.12 |
| **IBS** | +0.15 | +0.98 | -1.35 |
| **TSI** | +0.12 | -1.35 | +1.04 |

![spectrum](spectrum_fit.png)
