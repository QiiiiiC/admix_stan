# `(((IBS,TSI.1),EAS),TSI.2)`

**Normal, shared Ne, recent grid** | topology 19 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1112.16 | +- 0.45 (MC) |
| logZ (importance sampling) | -1081.84 | |
| ESS of the IS weights | 1.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 7 | 23 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 200.1 +- 0.7 | 211.1 |
| 3 | MERGE | n1 + EAS -> n2 | 208.2 +- 0.6 | 419.3 |
| 4 | MERGE | TSI.1 + n2 -> root | 486.4 +- 0.6 | 905.7 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `TSI.1`; 1.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 42,474,554 | 21,000,168 |
| `IBS` | 2,477,171 | 1,717,797 |
| `TSI` | 681,162 | 642,604 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 720,090 | 0.02 |
| `IBS` | 359,553 | 0.06 |
| `TSI` | 850,327 | 0.03 |
| `TSI.1` | 1,427,993,289 | 0.05 |
| `TSI.2` | 467,467 | 0.03 |
| `n1` | 1,175 | 0.02 |
| `n2` | 0 | 0.02 |
| `root` | 14,263 | 0.00 |

log-Ne random-walk step scale tau = 2.192

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -881.8 | 222 | 47.73 |
| SNP | +22.0 | 6 | 7.22 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +1.65 | -2.35 | -0.92 |
| **IBS** | -2.35 | +2.62 | +1.98 |
| **TSI** | -0.92 | +1.98 | -0.12 |

![spectrum](spectrum_fit.png)
