# `(((EAS.1,IBS),TSI),EAS.2)`

**Poisson, shared Ne, recent grid** | topology 05 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5210.50 | +- 0.48 (MC) |
| logZ (importance sampling) | -5179.12 | |
| ESS of the IS weights | 1.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 7 | 21 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 225.7 +- 1.1 | 236.7 |
| 3 | MERGE | n1 + TSI -> n2 | 1.0 +- 0.0 | 237.7 |
| 4 | MERGE | EAS.1 + n2 -> root | 158.8 +- 1.9 | 396.5 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 131,723,401 | 76,933,927 |
| `IBS` | 6,775,386 | 4,078,595 |
| `TSI` | 1,878,182 | 1,215,787 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 5,112,091 | 0.04 |
| `IBS` | 297,323 | 0.05 |
| `TSI` | 398,476 | 0.02 |
| `EAS.1` | 819,867 | 0.04 |
| `EAS.2` | 7,278 | 0.02 |
| `n1` | 788 | 0.02 |
| `n2` | 883 | 0.02 |
| `root` | 25 | 0.01 |

log-Ne random-walk step scale tau = 2.839

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,955.2 | 222 | 88.49 |
| SNP | +33.8 | 6 | 3.28 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.11 | +0.65 | -0.88 |
| **IBS** | +0.65 | -1.22 | -0.03 |
| **TSI** | -0.88 | -0.03 | +1.70 |

![spectrum](spectrum_fit.png)
