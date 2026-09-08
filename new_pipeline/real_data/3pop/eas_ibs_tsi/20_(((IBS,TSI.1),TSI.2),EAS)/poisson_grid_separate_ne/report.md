# `(((IBS,TSI.1),TSI.2),EAS)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 20 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -4965.93 | +- 0.35 (MC) |
| logZ (importance sampling) | -4932.90 | |
| ESS of the IS weights | 6.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 2 / 6 | 25 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 68.5 +- 4.7 | 68.5 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 6.5 +- 0.2 | 74.9 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 77.4 +- 11.1 | 152.3 |
| 4 | MERGE | n2 + EAS -> root | 151.3 +- 5.6 | 303.6 |

## Admixture fraction

**f = 0.988 +- 0.001** (fraction from `TSI.1`; 0.012 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 76,013,526 | 38,225,902 |
| `IBS` | 1,739,238 | 1,514,405 |
| `TSI` | 1,231,936 | 986,607 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 831,677 | 0.02 |
| `IBS` | 818,355 | 0.08 |
| `TSI` | 386,394 | 0.05 |
| `TSI.1` | 372,746 | 0.08 |
| `TSI.2` | 85,651 | 0.16 |
| `n1` | 69,957 | 0.15 |
| `n2` | 23,167 | 0.23 |
| `root` | 1,528 | 0.04 |

log-Ne random-walk step scale tau_ibd = 2.154

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,563 | 1,650 |
| `IBS` | 149,349 | 138,202 |
| `TSI` | 108,044 | 122,769 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,719 | 0.01 |
| `IBS` | 141,290 | 0.04 |
| `TSI` | 133,636 | 0.03 |
| `TSI.1` | 100,325 | 0.04 |
| `TSI.2` | 79,939 | 0.06 |
| `n1` | 82,725 | 0.06 |
| `n2` | 47,434 | 0.02 |
| `root` | 26,084 | 0.03 |

log-Ne random-walk step scale tau_snp = 0.635

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,711.8 | 222 | 38.55 |
| SNP | +35.0 | 6 | 2.86 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.07 | +0.05 | -0.19 |
| **IBS** | +0.05 | -0.26 | +0.17 |
| **TSI** | -0.19 | +0.17 | +0.21 |

![spectrum](spectrum_fit.png)
