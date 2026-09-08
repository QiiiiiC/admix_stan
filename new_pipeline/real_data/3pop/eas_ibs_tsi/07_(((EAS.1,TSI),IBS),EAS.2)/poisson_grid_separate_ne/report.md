# `(((EAS.1,TSI),IBS),EAS.2)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 07 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -815.24 | +- 0.26 (MC) |
| logZ (importance sampling) | -788.29 | |
| ESS of the IS weights | 1.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 6 | 26 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 130.5 +- 0.7 | 130.5 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 29.2 +- 1.0 | 159.6 |
| 3 | MERGE | n1 + IBS -> n2 | 4.1 +- 0.1 | 163.7 |
| 4 | MERGE | EAS.1 + n2 -> root | 118.1 +- 0.8 | 281.9 |

## Admixture fraction

**f = 0.998 +- 0.000** (fraction from `EAS.1`; 0.002 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,849,726 | 3,965,270 |
| `IBS` | 3,159,133 | 2,093,088 |
| `TSI` | 1,019,327 | 776,042 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,406,817 | 0.03 |
| `IBS` | 293,557 | 0.04 |
| `TSI` | 415,444 | 0.04 |
| `EAS.1` | 41,474 | 0.02 |
| `EAS.2` | 23,955 | 0.04 |
| `n1` | 24,877 | 0.03 |
| `n2` | 12,141 | 0.04 |
| `root` | 14,939 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.874

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 901 | 960 |
| `IBS` | 121,995 | 120,821 |
| `TSI` | 113,513 | 113,258 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,092 | 0.02 |
| `IBS` | 120,920 | 0.04 |
| `TSI` | 113,912 | 0.03 |
| `EAS.1` | 2,640 | 0.01 |
| `EAS.2` | 55,855 | 0.03 |
| `n1` | 56,983 | 0.03 |
| `n2` | 53,635 | 0.03 |
| `root` | 15,542 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.818

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -605.4 | 222 | 3.44 |
| SNP | +36.0 | 6 | 2.54 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.03 | +0.04 | +0.03 |
| **IBS** | +0.04 | -0.20 | +0.14 |
| **TSI** | +0.03 | +0.14 | -0.18 |

![spectrum](spectrum_fit.png)
