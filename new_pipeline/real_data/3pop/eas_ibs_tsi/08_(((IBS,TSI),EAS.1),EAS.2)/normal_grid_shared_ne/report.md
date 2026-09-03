# `(((IBS,TSI),EAS.1),EAS.2)`

**Normal, shared Ne, recent grid** | topology 08 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +3906.71 | +- 0.15 (MC) |
| logZ (importance sampling) | +3919.52 | |
| ESS of the IS weights | 2.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 13 | 24 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 131.9 +- 0.4 | 131.9 |
| 2 | MERGE | IBS + TSI -> n1 | 17.1 +- 0.7 | 148.9 |
| 3 | MERGE | EAS.2 + n1 -> n2 | 155.8 +- 2.0 | 304.8 |
| 4 | MERGE | EAS.1 + n2 -> root | 255.5 +- 0.8 | 560.3 |

## Admixture fraction

**f = 0.806 +- 0.001** (fraction from `EAS.1`; 0.194 from `EAS.2`)

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 4,627,924 | 4,072,906 |
| `IBS` | 2,158,099 | 1,478,195 |
| `TSI` | 894,165 | 712,488 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,398,713 | 0.03 |
| `IBS` | 274,608 | 0.03 |
| `TSI` | 406,220 | 0.02 |
| `EAS.1` | 36,634 | 0.02 |
| `EAS.2` | 5,025 | 0.04 |
| `n1` | 21,288 | 0.02 |
| `n2` | 1,004 | 0.01 |
| `root` | 15,005 | 0.00 |

log-Ne random-walk step scale tau = 1.843

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +4,033.5 | 222 | 3.13 |
| SNP | +30.8 | 6 | 4.29 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.37 | +0.50 | +0.23 |
| **IBS** | +0.50 | +1.92 | -3.08 |
| **TSI** | +0.23 | -3.08 | +2.47 |

![spectrum](spectrum_fit.png)
