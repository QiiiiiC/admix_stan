# `(((EAS,TSI.1),IBS),TSI.2)`

**Poisson, shared Ne, recent grid** | topology 17 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -4789.44 | +- 0.33 (MC) |
| logZ (importance sampling) | -4768.25 | |
| ESS of the IS weights | 1.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 1 / 6 | 32 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 84.8 +- 8.0 | 84.8 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 207.1 +- 8.2 | 291.9 |
| 3 | MERGE | n1 + IBS -> n2 | 10.6 +- 0.3 | 302.5 |
| 4 | MERGE | TSI.1 + n2 -> root | 1.0 +- 0.0 | 303.6 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 25,518,897 | 15,340,724 |
| `IBS` | 4,307,025 | 2,560,609 |
| `TSI` | 1,414,697 | 1,022,328 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,585,116 | 0.04 |
| `IBS` | 225,157 | 0.03 |
| `TSI` | 410,672 | 0.07 |
| `TSI.1` | 134,687 | 0.11 |
| `TSI.2` | 57 | 0.04 |
| `n1` | 59 | 0.03 |
| `n2` | 1,429 | 0.03 |
| `root` | 784 | 0.01 |

log-Ne random-walk step scale tau = 2.347

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,581.6 | 222 | 23093.30 |
| SNP | +32.8 | 6 | 3.62 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.40 | -0.88 | +0.10 |
| **IBS** | -0.88 | -0.28 | +2.09 |
| **TSI** | +0.10 | +2.09 | -2.15 |

![spectrum](spectrum_fit.png)
