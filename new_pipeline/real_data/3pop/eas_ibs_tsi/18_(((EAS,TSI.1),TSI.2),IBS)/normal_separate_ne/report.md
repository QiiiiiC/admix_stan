# `(((EAS,TSI.1),TSI.2),IBS)`

**Normal, separate IBD/SNP Ne** | topology 18 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +85.13 | +- 0.16 (MC) |
| logZ (importance sampling) | +103.01 | |
| ESS of the IS weights | 1.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 1 / 7 | 25 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 102.8 +- 0.3 | 102.8 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 26.6 +- 0.1 | 129.4 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 16.7 +- 0.1 | 146.1 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 147.2 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,671,389 | 0.03 |
| `IBS` | 219,026 | 0.02 |
| `TSI` | 430,731 | 0.03 |
| `TSI.1` | 76,886 | 0.02 |
| `TSI.2` | 22,299 | 0.02 |
| `n1` | 25,822 | 0.01 |
| `n2` | 119,380 | 0.02 |
| `root` | 75,987 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.542

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 725 | 0.01 |
| `IBS` | 116,294 | 0.14 |
| `TSI` | 134,572 | 0.07 |
| `TSI.1` | 54,236 | 0.05 |
| `TSI.2` | 10,246 | 0.01 |
| `n1` | 9,326 | 0.01 |
| `n2` | 18,536 | 0.01 |
| `root` | 18,956 | 0.01 |

log-Ne random-walk step scale tau_snp = 1.042

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +244.3 | 222 | 36.04 |
| SNP | +40.2 | 6 | 1.13 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.04 | -0.08 | +0.01 |
| **IBS** | -0.08 | -0.16 | +0.33 |
| **TSI** | +0.01 | +0.33 | -0.33 |

![spectrum](spectrum_fit.png)
