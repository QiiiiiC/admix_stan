# `(((IBS,TSI),EAS.1),EAS.2)`

**Normal, shared Ne, recent grid** | topology 08 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +3906.57 | +- 0.16 (MC) |
| logZ (importance sampling) | +3918.69 | |
| ESS of the IS weights | 6.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 1 / 0 | 116 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 131.8 +- 0.4 | 131.8 |
| 2 | MERGE | IBS + TSI -> n1 | 17.1 +- 0.7 | 148.9 |
| 3 | MERGE | EAS.2 + n1 -> n2 | 155.9 +- 2.2 | 304.8 |
| 4 | MERGE | EAS.1 + n2 -> root | 255.5 +- 0.8 | 560.4 |

## Admixture fraction

**f = 0.806 +- 0.001** (fraction from `EAS.1`; 0.194 from `EAS.2`)

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,630,203 | 4,074,648 |
| `IBS` | 2,160,023 | 1,479,516 |
| `TSI` | 894,001 | 712,503 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,400,080 | 0.03 |
| `IBS` | 274,593 | 0.03 |
| `TSI` | 406,564 | 0.02 |
| `EAS.1` | 36,675 | 0.02 |
| `EAS.2` | 5,023 | 0.04 |
| `n1` | 21,307 | 0.02 |
| `n2` | 1,005 | 0.01 |
| `root` | 15,005 | 0.00 |

log-Ne random-walk step scale tau = 1.843

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +4,033.2 | 222 | 3.13 |
| SNP | +30.9 | 6 | 4.25 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.35 | +0.48 | +0.20 |
| **IBS** | +0.48 | +1.95 | -3.06 |
| **TSI** | +0.20 | -3.06 | +2.49 |

![spectrum](spectrum_fit.png)
