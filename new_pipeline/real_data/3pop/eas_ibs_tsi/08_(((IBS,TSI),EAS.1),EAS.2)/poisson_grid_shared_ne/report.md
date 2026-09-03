# `(((IBS,TSI),EAS.1),EAS.2)`

**Poisson, shared Ne, recent grid** | topology 08 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5161.16 | +- 1.05 (MC) |
| logZ (importance sampling) | -5099.64 | |
| ESS of the IS weights | 3.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 7 | 21 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | IBS + TSI -> n1 | 156.1 +- 1.2 | 167.1 |
| 3 | MERGE | EAS.2 + n1 -> n2 | 203.3 +- 1.1 | 370.4 |
| 4 | MERGE | EAS.1 + n2 -> root | 31.6 +- 0.3 | 402.0 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 239,435,643 | 141,983,395 |
| `IBS` | 7,863,342 | 4,613,083 |
| `TSI` | 2,176,024 | 1,356,687 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 7,379,990 | 0.05 |
| `IBS` | 292,420 | 0.09 |
| `TSI` | 400,555 | 0.09 |
| `EAS.1` | 823,004 | 0.05 |
| `EAS.2` | 2,591,308 | 0.04 |
| `n1` | 11,488 | 0.03 |
| `n2` | 194 | 0.01 |
| `root` | 13 | 0.01 |

log-Ne random-walk step scale tau = 3.223

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,886.0 | 222 | 40.67 |
| SNP | +25.6 | 6 | 6.01 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.33 | +0.45 | +0.19 |
| **IBS** | +0.45 | +1.87 | -2.93 |
| **TSI** | +0.19 | -2.93 | +2.39 |

![spectrum](spectrum_fit.png)
