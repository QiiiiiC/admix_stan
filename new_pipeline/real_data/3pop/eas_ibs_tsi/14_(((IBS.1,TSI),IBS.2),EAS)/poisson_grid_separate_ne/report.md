# `(((IBS.1,TSI),IBS.2),EAS)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 14 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -4962.77 | +- 0.46 (MC) |
| logZ (importance sampling) | -4928.24 | |
| ESS of the IS weights | 1.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 3 | 24 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 77.9 +- 0.4 | 77.9 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 56.2 +- 0.7 | 134.0 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 90.9 +- 0.3 | 225.0 |
| 4 | MERGE | n2 + EAS -> root | 110.1 +- 0.8 | 335.1 |

## Admixture fraction

**f = 0.051 +- 0.001** (fraction from `IBS.1`; 0.949 from `IBS.2`)

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 62,766,867 | 35,030,655 |
| `IBS` | 1,958,059 | 1,637,712 |
| `TSI` | 1,914,608 | 1,249,309 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 833,155 | 0.03 |
| `IBS` | 799,175 | 0.06 |
| `TSI` | 390,290 | 0.10 |
| `IBS.1` | 18,689 | 0.05 |
| `IBS.2` | 53,772 | 0.03 |
| `n1` | 48,853 | 0.02 |
| `n2` | 6,593 | 0.02 |
| `root` | 439 | 0.01 |

log-Ne random-walk step scale tau_ibd = 2.297

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 2,055 | 1,971 |
| `IBS` | 95,603 | 90,788 |
| `TSI` | 95,368 | 99,048 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,899 | 0.02 |
| `IBS` | 86,974 | 0.02 |
| `TSI` | 95,882 | 0.03 |
| `IBS.1` | 31,438 | 0.02 |
| `IBS.2` | 89,829 | 0.02 |
| `n1` | 73,405 | 0.02 |
| `n2` | 37,681 | 0.02 |
| `root` | 15,956 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.724

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,710.3 | 222 | 38.51 |
| SNP | +35.2 | 6 | 2.80 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.03 | +0.00 | -0.07 |
| **IBS** | +0.00 | -0.24 | +0.26 |
| **TSI** | -0.07 | +0.26 | -0.11 |

![spectrum](spectrum_fit.png)
