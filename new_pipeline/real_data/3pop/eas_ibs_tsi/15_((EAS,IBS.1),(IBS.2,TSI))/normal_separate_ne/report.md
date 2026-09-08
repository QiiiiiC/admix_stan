# `((EAS,IBS.1),(IBS.2,TSI))`

**Normal, separate IBD/SNP Ne** | topology 15 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +4112.25 | +- 0.17 (MC) |
| logZ (importance sampling) | +4131.02 | |
| ESS of the IS weights | 1.3 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 1 / 1 | 26 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 71.5 +- 0.6 | 71.5 |
| 2 | MERGE | IBS.1 + EAS -> n1 | 58.1 +- 0.4 | 129.6 |
| 3 | MERGE | IBS.2 + TSI -> n2 | 11.6 +- 0.2 | 141.2 |
| 4 | MERGE | n1 + n2 -> root | 135.3 +- 1.9 | 276.6 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `IBS.1`; 1.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,544,482 | 0.04 |
| `IBS` | 916,293 | 0.05 |
| `TSI` | 404,620 | 0.05 |
| `IBS.1` | 42,484 | 0.02 |
| `IBS.2` | 73,706 | 0.03 |
| `n1` | 42,371 | 0.02 |
| `n2` | 29,502 | 0.03 |
| `root` | 15,408 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.399

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,051 | 0.01 |
| `IBS` | 118,447 | 0.04 |
| `TSI` | 95,302 | 0.04 |
| `IBS.1` | 2,742 | 0.01 |
| `IBS.2` | 95,054 | 0.04 |
| `n1` | 2,687 | 0.01 |
| `n2` | 59,733 | 0.04 |
| `root` | 14,101 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.858

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +4,251.0 | 222 | 1.03 |
| SNP | +39.4 | 6 | 1.41 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.07 | -0.08 | -0.06 |
| **IBS** | -0.08 | -0.07 | +0.23 |
| **TSI** | -0.06 | +0.23 | -0.10 |

![spectrum](spectrum_fit.png)
