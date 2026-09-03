# `(((IBS.1,TSI),IBS.2),EAS)`

**Poisson, shared Ne, recent grid** | topology 14 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5353.37 | +- 0.93 (MC) |
| logZ (importance sampling) | -5298.49 | |
| ESS of the IS weights | 1.4 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 13 | 21 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 243.4 +- 0.9 | 243.4 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 1.0 +- 0.0 | 244.4 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 1.0 +- 0.0 | 245.4 |
| 4 | MERGE | n2 + EAS -> root | 120.6 +- 2.6 | 366.0 |

## Admixture fraction

**f = 0.002 +- 0.000** (fraction from `IBS.1`; 0.998 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 80,743,378 | 38,572,738 |
| `IBS` | 4,714,357 | 2,903,937 |
| `TSI` | 1,486,877 | 1,034,221 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 829,764 | 0.03 |
| `IBS` | 299,257 | 0.12 |
| `TSI` | 401,418 | 0.08 |
| `IBS.1` | 682 | 0.03 |
| `IBS.2` | 811 | 0.05 |
| `n1` | 799 | 0.04 |
| `n2` | 676 | 0.03 |
| `root` | 131 | 0.08 |

log-Ne random-walk step scale tau = 2.445

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,091.4 | 222 | 112.02 |
| SNP | +31.5 | 6 | 4.04 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.10 | +0.65 | -0.85 |
| **IBS** | +0.65 | -1.20 | -0.04 |
| **TSI** | -0.85 | -0.04 | +1.65 |

![spectrum](spectrum_fit.png)
