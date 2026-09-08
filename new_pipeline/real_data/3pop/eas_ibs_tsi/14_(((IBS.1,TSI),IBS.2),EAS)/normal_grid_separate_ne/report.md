# `(((IBS.1,TSI),IBS.2),EAS)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 14 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +287.91 | +- 0.40 (MC) |
| logZ (importance sampling) | +336.03 | |
| ESS of the IS weights | 1.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 3 | 25 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 76.3 +- 0.9 | 76.3 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 60.0 +- 2.0 | 136.2 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 53.4 +- 2.5 | 189.6 |
| 4 | MERGE | n2 + EAS -> root | 119.7 +- 14.5 | 309.3 |

## Admixture fraction

**f = 0.033 +- 0.003** (fraction from `IBS.1`; 0.967 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 57,975,725 | 20,477,406 |
| `IBS` | 1,874,140 | 1,112,537 |
| `TSI` | 660,002 | 685,666 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 864,985 | 0.03 |
| `IBS` | 824,373 | 0.08 |
| `TSI` | 396,270 | 0.06 |
| `IBS.1` | 11,091 | 0.34 |
| `IBS.2` | 58,656 | 0.05 |
| `n1` | 43,595 | 0.08 |
| `n2` | 32,929 | 0.20 |
| `root` | 399 | 0.69 |

log-Ne random-walk step scale tau_ibd = 2.645

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,482 | 1,540 |
| `IBS` | 106,656 | 113,407 |
| `TSI` | 118,983 | 128,429 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,756 | 0.05 |
| `IBS` | 96,453 | 0.07 |
| `TSI` | 138,207 | 0.14 |
| `IBS.1` | 42,344 | 0.05 |
| `IBS.2` | 73,058 | 0.07 |
| `n1` | 58,808 | 0.05 |
| `n2` | 49,807 | 0.04 |
| `root` | 10,710 | 0.07 |

log-Ne random-walk step scale tau_snp = 0.809

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +542.2 | 222 | 34.34 |
| SNP | +37.4 | 6 | 2.06 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.03 | +0.27 | -0.32 |
| **IBS** | +0.27 | -0.34 | -0.18 |
| **TSI** | -0.32 | -0.18 | +0.79 |

![spectrum](spectrum_fit.png)
