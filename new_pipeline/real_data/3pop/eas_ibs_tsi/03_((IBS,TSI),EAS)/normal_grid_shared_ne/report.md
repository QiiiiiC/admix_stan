# `((IBS,TSI),EAS)`

**Normal, shared Ne, recent grid** | topology 03 of 21 | tree (no admixture) | 2 events, 5 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1104.61 | +- 0.52 (MC) |
| logZ (importance sampling) | -1074.97 | |
| ESS of the IS weights | 4.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -17.65 | already applied |
| mode kept / MAP start / runtime | 1 / 8 | 22 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | IBS + TSI -> n1 | 211.3 +- 0.5 | 211.3 |
| 2 | MERGE | EAS + n1 -> root | 207.0 +- 1.2 | 418.4 |

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 97,482,890 | 44,741,168 |
| `IBS` | 2,094,893 | 1,528,605 |
| `TSI` | 516,795 | 652,264 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 718,249 | 0.03 |
| `IBS` | 359,543 | 0.03 |
| `TSI` | 470,103 | 0.05 |
| `n1` | 1,162 | 0.02 |
| `root` | 3 | 0.02 |

log-Ne random-walk step scale tau = 2.120

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -897.9 | 222 | 47.87 |
| SNP | +14.4 | 6 | 9.73 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +1.30 | -1.16 | -1.41 |
| **IBS** | -1.16 | +3.32 | -1.20 |
| **TSI** | -1.41 | -1.20 | +3.81 |

![spectrum](spectrum_fit.png)
