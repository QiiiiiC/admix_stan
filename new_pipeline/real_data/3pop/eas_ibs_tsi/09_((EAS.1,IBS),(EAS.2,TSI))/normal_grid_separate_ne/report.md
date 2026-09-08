# `((EAS.1,IBS),(EAS.2,TSI))`

**Normal, separate IBD/SNP Ne, recent grid** | topology 09 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +451.20 | +- 0.64 (MC) |
| logZ (importance sampling) | +497.86 | |
| ESS of the IS weights | 1.9 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 11 | 27 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 124.5 +- 0.6 | 124.5 |
| 2 | MERGE | EAS.1 + IBS -> n1 | 5.7 +- 0.3 | 130.2 |
| 3 | MERGE | EAS.2 + TSI -> n2 | 16.8 +- 0.3 | 147.0 |
| 4 | MERGE | n1 + n2 -> root | 1.0 +- 0.0 | 148.0 |

## Admixture fraction

**f = 0.001 +- 0.000** (fraction from `EAS.1`; 0.999 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 3,698,194 | 3,445,901 |
| `IBS` | 5,345,213 | 3,871,397 |
| `TSI` | 3,402,828 | 2,112,157 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,637,127 | 0.06 |
| `IBS` | 555,245 | 0.08 |
| `TSI` | 316,487 | 0.07 |
| `EAS.1` | 6,335 | 0.05 |
| `EAS.2` | 37,305 | 0.04 |
| `n1` | 6,398 | 0.05 |
| `n2` | 225,306 | 0.02 |
| `root` | 74,072 | 0.03 |

log-Ne random-walk step scale tau_ibd = 2.290

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 707 | 714 |
| `IBS` | 227,486 | 222,170 |
| `TSI` | 144,975 | 142,460 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 694 | 0.02 |
| `IBS` | 231,771 | 0.14 |
| `TSI` | 137,340 | 0.16 |
| `EAS.1` | 23,111 | 0.10 |
| `EAS.2` | 12,937 | 0.05 |
| `n1` | 23,703 | 0.10 |
| `n2` | 18,851 | 0.07 |
| `root` | 19,126 | 0.07 |

log-Ne random-walk step scale tau_snp = 1.564

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +722.4 | 222 | 32.01 |
| SNP | +31.8 | 6 | 3.95 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.17 | +0.37 | -0.04 |
| **IBS** | +0.37 | -0.23 | -0.51 |
| **TSI** | -0.04 | -0.51 | +0.55 |

![spectrum](spectrum_fit.png)
