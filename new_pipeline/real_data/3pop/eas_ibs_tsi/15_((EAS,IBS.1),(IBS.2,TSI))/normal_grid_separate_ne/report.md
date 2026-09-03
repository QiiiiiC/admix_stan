# `((EAS,IBS.1),(IBS.2,TSI))`

**Normal, separate IBD/SNP Ne, recent grid** | topology 15 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +823.64 | +- 0.85 (MC) |
| logZ (importance sampling) | +877.34 | |
| ESS of the IS weights | 1.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 7 | 30 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 34.7 +- 0.5 | 34.7 |
| 2 | MERGE | IBS.1 + EAS -> n1 | 30.0 +- 0.8 | 64.7 |
| 3 | MERGE | IBS.2 + TSI -> n2 | 1.1 +- 0.0 | 65.8 |
| 4 | MERGE | n1 + n2 -> root | 374.1 +- 1.5 | 439.9 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `IBS.1`; 1.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 6,135,366 | 5,907,354 |
| `IBS` | 571,687 | 488,731 |
| `TSI` | 712,969 | 611,363 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 4,012,485 | 0.07 |
| `IBS` | 306,823 | 0.04 |
| `TSI` | 429,176 | 0.04 |
| `IBS.1` | 261,641 | 0.05 |
| `IBS.2` | 168,132 | 0.04 |
| `n1` | 261,633 | 0.05 |
| `n2` | 217,342 | 0.03 |
| `root` | 2 | 0.07 |

log-Ne random-walk step scale tau_ibd = 0.896

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 5,027 | 5,028 |
| `IBS` | 5,205 | 5,204 |
| `TSI` | 5,200 | 5,199 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 5,032 | 0.01 |
| `IBS` | 5,203 | 0.01 |
| `TSI` | 5,199 | 0.01 |
| `IBS.1` | 5,053 | 0.01 |
| `IBS.2` | 5,202 | 0.01 |
| `n1` | 5,053 | 0.01 |
| `n2` | 5,193 | 0.01 |
| `root` | 5,164 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.004

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +2,668.1 | 222 | 14.21 |
| SNP | -1,518.3 | 6 | 520.63 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +10.12 | -9.91 | -10.10 |
| **IBS** | -9.91 | -21.60 | +43.38 |
| **TSI** | -10.10 | +43.38 | -21.56 |

![spectrum](spectrum_fit.png)
