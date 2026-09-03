# `(((EAS.1,IBS),TSI),EAS.2)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 05 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5169.02 | +- 0.29 (MC) |
| logZ (importance sampling) | -5144.89 | |
| ESS of the IS weights | 2.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 7 | 24 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 152.1 +- 1.9 | 163.2 |
| 3 | MERGE | n1 + TSI -> n2 | 1.0 +- 0.0 | 164.2 |
| 4 | MERGE | EAS.1 + n2 -> root | 141.4 +- 1.4 | 305.5 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 120,095,286 | 77,770,534 |
| `IBS` | 5,574,456 | 3,553,379 |
| `TSI` | 1,779,195 | 1,164,240 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 5,655,937 | 0.02 |
| `IBS` | 299,209 | 0.05 |
| `TSI` | 402,509 | 0.02 |
| `EAS.1` | 819,457 | 0.01 |
| `EAS.2` | 12,378 | 0.02 |
| `n1` | 6,517 | 0.02 |
| `n2` | 12,458 | 0.03 |
| `root` | 1,482 | 0.01 |

log-Ne random-walk step scale tau_ibd = 2.906

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 2,642 | 2,522 |
| `IBS` | 225,793 | 215,929 |
| `TSI` | 113,873 | 109,596 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,992 | 0.01 |
| `IBS` | 171,678 | 0.05 |
| `TSI` | 99,358 | 0.03 |
| `EAS.1` | 1,715 | 0.01 |
| `EAS.2` | 37,779 | 0.02 |
| `n1` | 32,216 | 0.01 |
| `n2` | 32,184 | 0.01 |
| `root` | 11,401 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.775

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,846.4 | 222 | 40.14 |
| SNP | +40.5 | 6 | 1.02 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.02 | -0.24 | +0.28 |
| **IBS** | -0.24 | +0.61 | -0.17 |
| **TSI** | +0.28 | -0.17 | -0.38 |

![spectrum](spectrum_fit.png)
