# `(((EAS.1,TSI),IBS),EAS.2)`

**Poisson, shared Ne, recent grid** | topology 07 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5272.71 | +- 0.91 (MC) |
| logZ (importance sampling) | -5220.96 | |
| ESS of the IS weights | 1.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 1 | 20 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 236.0 +- 0.7 | 247.0 |
| 3 | MERGE | n1 + IBS -> n2 | 1.0 +- 0.0 | 248.0 |
| 4 | MERGE | EAS.1 + n2 -> root | 106.9 +- 1.5 | 354.9 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 80,281,976 | 40,560,817 |
| `IBS` | 2,408,988 | 1,596,414 |
| `TSI` | 779,083 | 651,143 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,146,469 | 0.04 |
| `IBS` | 293,945 | 0.03 |
| `TSI` | 415,712 | 0.19 |
| `EAS.1` | 821,472 | 0.02 |
| `EAS.2` | 794 | 0.02 |
| `n1` | 1,009 | 0.02 |
| `n2` | 593 | 0.01 |
| `root` | 216 | 0.03 |

log-Ne random-walk step scale tau = 2.690

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,051.9 | 222 | 127.46 |
| SNP | +39.2 | 6 | 1.49 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.02 | -0.33 | +0.30 |
| **IBS** | -0.33 | +1.02 | -0.42 |
| **TSI** | +0.30 | -0.42 | -0.18 |

![spectrum](spectrum_fit.png)
