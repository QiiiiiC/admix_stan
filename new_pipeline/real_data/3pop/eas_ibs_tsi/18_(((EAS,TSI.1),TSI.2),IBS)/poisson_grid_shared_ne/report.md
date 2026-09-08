# `(((EAS,TSI.1),TSI.2),IBS)`

**Poisson, shared Ne, recent grid** | topology 18 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -4815.81 | +- 0.29 (MC) |
| logZ (importance sampling) | -4796.41 | |
| ESS of the IS weights | 3.4 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 3 / 7 | 29 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 286.0 +- 0.4 | 286.0 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 6.5 +- 0.1 | 292.5 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 9.4 +- 0.1 | 301.9 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 302.9 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 24,735,760 | 18,189,388 |
| `IBS` | 5,891,129 | 3,393,583 |
| `TSI` | 3,270,234 | 2,019,873 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,581,065 | 0.03 |
| `IBS` | 225,040 | 0.03 |
| `TSI` | 308,034 | 0.07 |
| `TSI.1` | 14,757 | 0.03 |
| `TSI.2` | 46 | 0.01 |
| `n1` | 52 | 0.01 |
| `n2` | 1,860 | 0.02 |
| `root` | 798 | 0.02 |

log-Ne random-walk step scale tau = 2.728

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,614.0 | 222 | 22005.95 |
| SNP | +36.6 | 6 | 2.35 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.45 | -0.41 | -0.49 |
| **IBS** | -0.41 | -1.21 | +2.12 |
| **TSI** | -0.49 | +2.12 | -1.06 |

![spectrum](spectrum_fit.png)
