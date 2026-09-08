# `(((EAS,TSI.1),IBS),TSI.2)`

**Poisson, separate IBD/SNP Ne** | topology 17 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -4245.15 | +- 0.16 (MC) |
| logZ (importance sampling) | -4228.04 | |
| ESS of the IS weights | 2.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 1 / 7 | 19 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 109.0 +- 0.5 | 109.0 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 21.3 +- 0.1 | 130.2 |
| 3 | MERGE | n1 + IBS -> n2 | 50.6 +- 0.9 | 180.8 |
| 4 | MERGE | TSI.1 + n2 -> root | 1.1 +- 0.0 | 181.9 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,577,258 | 0.02 |
| `IBS` | 235,071 | 0.01 |
| `TSI` | 440,896 | 0.06 |
| `TSI.1` | 62,439 | 0.02 |
| `TSI.2` | 36,837 | 0.02 |
| `n1` | 37,749 | 0.02 |
| `n2` | 83,970 | 0.01 |
| `root` | 65,245 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.365

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 762 | 0.01 |
| `IBS` | 128,239 | 0.06 |
| `TSI` | 179,853 | 0.09 |
| `TSI.1` | 79,762 | 0.07 |
| `TSI.2` | 4,927 | 0.01 |
| `n1` | 5,434 | 0.01 |
| `n2` | 18,689 | 0.01 |
| `root` | 19,148 | 0.01 |

log-Ne random-walk step scale tau_snp = 1.024

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,097.0 | 222 | 81.64 |
| SNP | +39.3 | 6 | 1.45 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.04 | -0.06 | -0.03 |
| **IBS** | -0.06 | -0.32 | +0.47 |
| **TSI** | -0.03 | +0.47 | -0.39 |

![spectrum](spectrum_fit.png)
