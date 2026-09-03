# `(((IBS,TSI.1),TSI.2),EAS)`

**Poisson, shared Ne, recent grid** | topology 20 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5347.04 | +- 0.77 (MC) |
| logZ (importance sampling) | -5299.11 | |
| ESS of the IS weights | 3.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 13 | 22 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 236.9 +- 0.6 | 247.9 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 1.0 +- 0.0 | 248.9 |
| 4 | MERGE | n2 + EAS -> root | 102.5 +- 1.3 | 351.4 |

## Admixture fraction

**f = 0.831 +- 0.003** (fraction from `TSI.1`; 0.169 from `TSI.2`)

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 79,529,945 | 40,893,904 |
| `IBS` | 4,743,863 | 2,992,978 |
| `TSI` | 1,591,857 | 1,471,161 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 828,617 | 0.03 |
| `IBS` | 296,324 | 0.05 |
| `TSI` | 1,232,469 | 0.27 |
| `TSI.1` | 2,506,664 | 0.30 |
| `TSI.2` | 12,643 | 0.04 |
| `n1` | 656 | 0.02 |
| `n2` | 570 | 0.03 |
| `root` | 250 | 0.01 |

log-Ne random-walk step scale tau = 2.544

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,084.6 | 222 | 130.53 |
| SNP | +24.2 | 6 | 6.46 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.04 | +0.82 | -0.91 |
| **IBS** | +0.82 | -1.28 | -0.30 |
| **TSI** | -0.91 | -0.30 | +2.03 |

![spectrum](spectrum_fit.png)
