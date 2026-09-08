# `(((EAS.1,IBS),TSI),EAS.2)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 05 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +4045.63 | +- 0.40 (MC) |
| logZ (importance sampling) | +4077.42 | |
| ESS of the IS weights | 1.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 7 | 27 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 126.6 +- 0.6 | 126.6 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 2.4 +- 0.1 | 129.0 |
| 3 | MERGE | n1 + TSI -> n2 | 16.2 +- 1.1 | 145.2 |
| 4 | MERGE | EAS.1 + n2 -> root | 142.6 +- 1.8 | 287.8 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 5,842,575 | 5,022,612 |
| `IBS` | 1,957,123 | 1,584,040 |
| `TSI` | 1,193,958 | 892,780 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,417,100 | 0.02 |
| `IBS` | 560,657 | 0.07 |
| `TSI` | 390,360 | 0.04 |
| `EAS.1` | 47,404 | 0.02 |
| `EAS.2` | 7,835 | 0.07 |
| `n1` | 8,008 | 0.07 |
| `n2` | 26,138 | 0.02 |
| `root` | 9,942 | 0.02 |

log-Ne random-walk step scale tau_ibd = 1.710

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,374 | 1,308 |
| `IBS` | 191,046 | 194,454 |
| `TSI` | 84,696 | 83,986 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,319 | 0.02 |
| `IBS` | 188,553 | 0.06 |
| `TSI` | 86,421 | 0.06 |
| `EAS.1` | 1,951 | 0.02 |
| `EAS.2` | 84,540 | 0.05 |
| `n1` | 82,257 | 0.05 |
| `n2` | 66,666 | 0.05 |
| `root` | 13,600 | 0.02 |

log-Ne random-walk step scale tau_snp = 0.843

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +4,247.1 | 222 | 1.17 |
| SNP | +33.8 | 6 | 3.27 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.17 | -0.16 | +0.49 |
| **IBS** | -0.16 | +0.66 | -0.37 |
| **TSI** | +0.49 | -0.37 | -0.59 |

![spectrum](spectrum_fit.png)
