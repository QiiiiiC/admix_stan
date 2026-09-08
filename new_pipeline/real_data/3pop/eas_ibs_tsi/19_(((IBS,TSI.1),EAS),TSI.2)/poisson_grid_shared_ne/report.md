# `(((IBS,TSI.1),EAS),TSI.2)`

**Poisson, shared Ne, recent grid** | topology 19 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5247.61 | +- 0.26 (MC) |
| logZ (importance sampling) | -5227.26 | |
| ESS of the IS weights | 1.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 1 / 9 | 31 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 85.6 +- 0.4 | 85.6 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 149.4 +- 0.3 | 235.0 |
| 3 | MERGE | n1 + EAS -> n2 | 170.2 +- 0.2 | 405.2 |
| 4 | MERGE | TSI.1 + n2 -> root | 221.9 +- 0.3 | 627.0 |

## Admixture fraction

**f = 0.009 +- 0.000** (fraction from `TSI.1`; 0.991 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 202,938,866 | 81,425,614 |
| `IBS` | 9,418,912 | 5,209,648 |
| `TSI` | 2,167,620 | 1,482,082 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 831,763 | 0.02 |
| `IBS` | 296,672 | 0.02 |
| `TSI` | 348,918 | 0.06 |
| `TSI.1` | 65,298 | 0.01 |
| `TSI.2` | 3,924,687,328 | 0.07 |
| `n1` | 940 | 0.01 |
| `n2` | 8 | 0.03 |
| `root` | 17,812 | 0.00 |

log-Ne random-walk step scale tau = 3.002

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,019.3 | 222 | 80.54 |
| SNP | +32.8 | 6 | 3.61 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.20 | +0.58 | -0.98 |
| **IBS** | +0.58 | -1.84 | +0.79 |
| **TSI** | -0.98 | +0.79 | +1.12 |

![spectrum](spectrum_fit.png)
