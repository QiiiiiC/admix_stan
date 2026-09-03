# `((EAS,TSI.1),(IBS,TSI.2))`

**Normal, separate IBD/SNP Ne, recent grid** | topology 21 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -2691.69 | +- 4.52 (MC) |
| logZ (importance sampling) | -2483.81 | |
| ESS of the IS weights | 1.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 13 | 24 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | TSI.1 + EAS -> n1 | 1.0 +- 0.0 | 12.0 |
| 3 | MERGE | TSI.2 + IBS -> n2 | 1.0 +- 0.0 | 13.0 |
| 4 | MERGE | n1 + n2 -> root | 338.4 +- 1.6 | 351.4 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `TSI.1`; 1.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 65,282,629 | 52,507,968 |
| `IBS` | 3,137,246 | 2,053,543 |
| `TSI` | 2,355,310 | 1,606,780 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 11,433,235 | 0.10 |
| `IBS` | 573,363 | 0.09 |
| `TSI` | 607,187 | 0.09 |
| `TSI.1` | 904,958 | 0.09 |
| `TSI.2` | 494,481 | 0.09 |
| `n1` | 904,975 | 0.09 |
| `n2` | 350,563 | 0.09 |
| `root` | 34 | 0.08 |

log-Ne random-walk step scale tau_ibd = 2.636

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 2,380 | 2,382 |
| `IBS` | 9,880 | 9,934 |
| `TSI` | 9,935 | 9,979 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,392 | 0.07 |
| `IBS` | 10,056 | 0.10 |
| `TSI` | 10,080 | 0.10 |
| `TSI.1` | 2,397 | 0.07 |
| `TSI.2` | 10,105 | 0.10 |
| `n1` | 2,397 | 0.07 |
| `n2` | 10,169 | 0.10 |
| `root` | 15,364 | 0.02 |

log-Ne random-walk step scale tau_snp = 0.643

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -2,264.2 | 222 | 56.92 |
| SNP | -26.2 | 6 | 23.28 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.39 | +0.45 | +0.32 |
| **IBS** | +0.45 | -0.39 | -0.50 |
| **TSI** | +0.32 | -0.50 | -0.14 |

![spectrum](spectrum_fit.png)
