# `(((EAS.1,TSI),EAS.2),IBS)`

**Poisson, shared Ne, recent grid** | topology 06 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -42160.82 | +- 0.57 (MC) |
| logZ (importance sampling) | -42126.97 | |
| ESS of the IS weights | 5.9 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 1 | 18 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 277.9 +- 1.3 | 288.9 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 1.1 +- 0.0 | 289.9 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 290.9 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 122,106,373 | 65,774,270 |
| `IBS` | 10,003,330 | 5,476,015 |
| `TSI` | 5,386,398 | 2,961,824 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 4,288,623 | 0.09 |
| `IBS` | 229,505 | 0.04 |
| `TSI` | 331,976 | 0.08 |
| `EAS.1` | 1,042,180 | 0.03 |
| `EAS.2` | 76,527 | 0.57 |
| `n1` | 53 | 0.05 |
| `n2` | 56 | 0.02 |
| `root` | 692 | 0.05 |

log-Ne random-walk step scale tau = 3.213

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,059.8 | 222 | 6110.24 |
| SNP | -34,818.5 | 6 | 11620.71 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +111.05 | -110.06 | -109.44 |
| **IBS** | -110.06 | +54.85 | +165.87 |
| **TSI** | -109.44 | +165.87 | +52.63 |

![spectrum](spectrum_fit.png)
