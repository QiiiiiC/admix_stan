# `(((EAS,TSI.1),IBS),TSI.2)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 17 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +54.96 | +- 0.50 (MC) |
| logZ (importance sampling) | +96.17 | |
| ESS of the IS weights | 1.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 7 | 25 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 132.0 +- 0.6 | 132.0 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 1.0 +- 0.0 | 133.0 |
| 3 | MERGE | n1 + IBS -> n2 | 12.4 +- 0.3 | 145.4 |
| 4 | MERGE | TSI.1 + n2 -> root | 1.1 +- 0.0 | 146.4 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,523,257 | 4,203,380 |
| `IBS` | 5,748,606 | 2,940,613 |
| `TSI` | 1,844,403 | 1,244,193 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,520,579 | 0.07 |
| `IBS` | 210,149 | 0.05 |
| `TSI` | 394,946 | 0.05 |
| `TSI.1` | 20,396 | 0.10 |
| `TSI.2` | 18,525 | 0.05 |
| `n1` | 18,521 | 0.05 |
| `n2` | 205,586 | 0.08 |
| `root` | 78,455 | 0.01 |

log-Ne random-walk step scale tau_ibd = 2.567

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 785 | 797 |
| `IBS` | 188,488 | 157,129 |
| `TSI` | 217,333 | 211,367 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 740 | 0.02 |
| `IBS` | 112,326 | 0.05 |
| `TSI` | 187,748 | 0.06 |
| `TSI.1` | 20,655 | 0.13 |
| `TSI.2` | 8,178 | 0.16 |
| `n1` | 8,177 | 0.16 |
| `n2` | 24,627 | 0.03 |
| `root` | 22,491 | 0.05 |

log-Ne random-walk step scale tau_snp = 1.089

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +344.7 | 222 | 35.27 |
| SNP | +33.3 | 6 | 3.43 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.01 | -0.02 | +0.04 |
| **IBS** | -0.02 | +0.03 | -0.00 |
| **TSI** | +0.04 | -0.00 | -0.07 |

![spectrum](spectrum_fit.png)
