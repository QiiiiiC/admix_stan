# `(((IBS,TSI),EAS.1),EAS.2)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 08 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -801.08 | +- 0.42 (MC) |
| logZ (importance sampling) | -775.85 | |
| ESS of the IS weights | 9.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 2 / 11 | 27 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 131.4 +- 1.4 | 131.4 |
| 2 | MERGE | IBS + TSI -> n1 | 31.3 +- 0.3 | 162.7 |
| 3 | MERGE | EAS.2 + n1 -> n2 | 95.8 +- 2.0 | 258.5 |
| 4 | MERGE | EAS.1 + n2 -> root | 17.3 +- 0.4 | 275.8 |

## Admixture fraction

**f = 0.890 +- 0.001** (fraction from `EAS.1`; 0.110 from `EAS.2`)

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,575,164 | 3,800,349 |
| `IBS` | 2,956,870 | 2,037,105 |
| `TSI` | 1,079,251 | 856,502 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,403,625 | 0.03 |
| `IBS` | 294,627 | 0.03 |
| `TSI` | 408,224 | 0.07 |
| `EAS.1` | 32,156 | 0.06 |
| `EAS.2` | 152,969 | 0.08 |
| `n1` | 12,371 | 0.04 |
| `n2` | 20,181 | 0.05 |
| `root` | 17,559 | 0.05 |

log-Ne random-walk step scale tau_ibd = 1.907

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 972 | 970 |
| `IBS` | 146,397 | 135,578 |
| `TSI` | 100,140 | 101,611 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,055 | 0.02 |
| `IBS` | 124,501 | 0.07 |
| `TSI` | 110,652 | 0.09 |
| `EAS.1` | 2,213 | 0.01 |
| `EAS.2` | 16,652 | 0.03 |
| `n1` | 46,688 | 0.05 |
| `n2` | 12,109 | 0.01 |
| `root` | 11,858 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.850

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -605.9 | 222 | 3.42 |
| SNP | +31.6 | 6 | 4.01 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.02 | -0.01 | +0.05 |
| **IBS** | -0.01 | -0.08 | +0.10 |
| **TSI** | +0.05 | +0.10 | -0.19 |

![spectrum](spectrum_fit.png)
