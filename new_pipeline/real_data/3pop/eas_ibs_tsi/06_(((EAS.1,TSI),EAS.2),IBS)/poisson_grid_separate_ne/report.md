# `(((EAS.1,TSI),EAS.2),IBS)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 06 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -7319.75 | +- 0.15 (MC) |
| logZ (importance sampling) | -7300.72 | |
| ESS of the IS weights | 1.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 7 | 22 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 226.9 +- 0.6 | 237.9 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 1.3 +- 0.0 | 239.2 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 240.2 |

## Admixture fraction

**f = 0.999 +- 0.000** (fraction from `EAS.1`; 0.001 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 110,840,430 | 76,916,686 |
| `IBS` | 6,769,168 | 4,172,645 |
| `TSI` | 3,194,017 | 2,079,866 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 6,594,123 | 0.02 |
| `IBS` | 230,289 | 0.02 |
| `TSI` | 333,153 | 0.06 |
| `EAS.1` | 1,019,803 | 0.02 |
| `EAS.2` | 1,513 | 0.06 |
| `n1` | 335 | 0.01 |
| `n2` | 625 | 0.01 |
| `root` | 4,287 | 0.01 |

log-Ne random-walk step scale tau_ibd = 3.260

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,748 | 1,685 |
| `IBS` | 198,879 | 191,697 |
| `TSI` | 234,920 | 226,494 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,432 | 0.01 |
| `IBS` | 158,233 | 0.07 |
| `TSI` | 198,335 | 0.05 |
| `EAS.1` | 1,309 | 0.01 |
| `EAS.2` | 30,023 | 0.03 |
| `n1` | 20,203 | 0.01 |
| `n2` | 20,050 | 0.01 |
| `root` | 22,651 | 0.01 |

log-Ne random-walk step scale tau_snp = 1.295

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -6,988.4 | 222 | 355.43 |
| SNP | +41.7 | 6 | 0.64 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.02 | +0.18 | -0.23 |
| **IBS** | +0.18 | -0.52 | +0.19 |
| **TSI** | -0.23 | +0.19 | +0.27 |

![spectrum](spectrum_fit.png)
