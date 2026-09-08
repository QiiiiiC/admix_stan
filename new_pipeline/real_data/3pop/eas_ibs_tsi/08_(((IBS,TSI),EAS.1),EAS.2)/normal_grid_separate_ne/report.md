# `(((IBS,TSI),EAS.1),EAS.2)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 08 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +3828.41 | +- 0.12 (MC) |
| logZ (importance sampling) | +3843.59 | |
| ESS of the IS weights | 1.9 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 11 | 27 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 132.1 +- 0.5 | 132.1 |
| 2 | MERGE | IBS + TSI -> n1 | 15.6 +- 0.3 | 147.7 |
| 3 | MERGE | EAS.2 + n1 -> n2 | 158.7 +- 2.2 | 306.5 |
| 4 | MERGE | EAS.1 + n2 -> root | 340.2 +- 1.0 | 646.6 |

## Admixture fraction

**f = 0.739 +- 0.001** (fraction from `EAS.1`; 0.261 from `EAS.2`)

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,747,523 | 4,004,909 |
| `IBS` | 1,835,529 | 1,536,306 |
| `TSI` | 866,637 | 731,040 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,397,194 | 0.03 |
| `IBS` | 274,567 | 0.03 |
| `TSI` | 406,760 | 0.03 |
| `EAS.1` | 179,233 | 0.02 |
| `EAS.2` | 3,015 | 0.03 |
| `n1` | 21,833 | 0.03 |
| `n2` | 1,263 | 0.03 |
| `root` | 18,354 | 0.00 |

log-Ne random-walk step scale tau_ibd = 1.815

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 3,239 | 3,369 |
| `IBS` | 100,552 | 102,985 |
| `TSI` | 99,043 | 99,837 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 3,387 | 0.01 |
| `IBS` | 107,810 | 0.04 |
| `TSI` | 100,013 | 0.04 |
| `EAS.1` | 2,052 | 0.01 |
| `EAS.2` | 44,761 | 0.03 |
| `n1` | 105,634 | 0.04 |
| `n2` | 74,465 | 0.03 |
| `root` | 14,194 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.629

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +4,033.5 | 222 | 3.12 |
| SNP | +41.8 | 6 | 0.60 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.01 | +0.00 | -0.02 |
| **IBS** | +0.00 | -0.21 | +0.22 |
| **TSI** | -0.02 | +0.22 | -0.16 |

![spectrum](spectrum_fit.png)
