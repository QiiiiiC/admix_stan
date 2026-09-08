# `(((EAS.1,TSI),IBS),EAS.2)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 07 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +3808.88 | +- 0.15 (MC) |
| logZ (importance sampling) | +3827.83 | |
| ESS of the IS weights | 1.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 1 / 6 | 27 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 127.7 +- 0.5 | 127.7 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 12.4 +- 0.3 | 140.2 |
| 3 | MERGE | n1 + IBS -> n2 | 12.9 +- 0.2 | 153.1 |
| 4 | MERGE | EAS.1 + n2 -> root | 130.8 +- 0.4 | 283.9 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,546,554 | 4,189,024 |
| `IBS` | 2,462,577 | 1,664,514 |
| `TSI` | 1,113,339 | 827,224 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,432,735 | 0.05 |
| `IBS` | 274,645 | 0.02 |
| `TSI` | 411,863 | 0.03 |
| `EAS.1` | 45,601 | 0.02 |
| `EAS.2` | 121,838 | 0.02 |
| `n1` | 109,770 | 0.02 |
| `n2` | 17,407 | 0.01 |
| `root` | 11,496 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.921

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,018 | 1,085 |
| `IBS` | 100,283 | 105,581 |
| `TSI` | 112,132 | 117,284 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,205 | 0.01 |
| `IBS` | 110,394 | 0.04 |
| `TSI` | 111,779 | 0.04 |
| `EAS.1` | 2,216 | 0.01 |
| `EAS.2` | 54,522 | 0.03 |
| `n1` | 61,088 | 0.03 |
| `n2` | 48,278 | 0.03 |
| `root` | 15,692 | 0.00 |

log-Ne random-walk step scale tau_snp = 0.876

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +4,030.7 | 222 | 3.12 |
| SNP | +40.9 | 6 | 0.92 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.00 | +0.03 | -0.03 |
| **IBS** | +0.03 | -0.26 | +0.21 |
| **TSI** | -0.03 | +0.21 | -0.15 |

![spectrum](spectrum_fit.png)
