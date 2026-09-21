# Growth + asymmetric-Ne b-loop: shared vs separate Ne

Overlapping 30-of-50 block subsamples; no independent-replicate error bars. Evidence is compared only
within a model variant and IBD source (graph scores average over event orders). Under `growth` every
candidate is misspecified in Ne, so times and sizes are descriptive, not recovery targets.

## Topology selection (ELBO ranking, search including the explicit-loop graph)

Ranked by ELBO averaged over event orders. The importance-sampling logZ is recorded but not used to rank:
pooled ESS is O(1–20) of 12,000 draws and Pareto k > 1 on essentially every fit, so it is one draw's weight.

| Scenario | IBD | Model | Complete | Reliable | Backbone correct | Explicit loop wins | Median rank of g15 | Median rank of g22 |
|---|---|---|---:|---:|---:|---:|---:|---:|
| growth | true | poisson_shared | 10 | 0 | 100% | 0% | 1 | 2 |
| growth | true | poisson_separate | 10 | 0 | 90% | 70% | 3 | 1 |
| growth | hapibd | poisson_shared | 10 | 0 | 100% | 0% | 1 | 2 |
| growth | hapibd | poisson_separate | 10 | 0 | 70% | 50% | 4 | 2 |
| constant | true | poisson_shared | 10 | 0 | 100% | 0% | 1 | 2 |
| constant | true | poisson_separate | 10 | 0 | 90% | 90% | 3 | 1 |
| constant | hapibd | poisson_shared | 10 | 0 | 100% | 0% | 1 | 2 |
| constant | hapibd | poisson_separate | 10 | 0 | 100% | 100% | 2 | 1 |

## Does a separate Ne help?  Paired on identical data

dELBO = separate − shared for the same graph on the same subsample; positive favours separate. Constants
are normalised so the extra parameters are paid for. 'Rescue' = separate picks a correct backbone where
shared did not; 'loss' the reverse. (dlogZ is in paired_ne_comparison.csv, subject to the caveat above.)

| Scenario | IBD | n | mean dELBO backbone | share > 0 | mean dELBO explicit | share > 0 | rescue | loss |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| growth | true | 10 | -57.2 | 0% | -32.9 | 0% | 0 | 1 |
| growth | hapibd | 10 | -69.0 | 0% | -37.4 | 0% | 0 | 3 |
| constant | true | 10 | -92.9 | 0% | -46.9 | 0% | 0 | 1 |
| constant | hapibd | 10 | -75.1 | 0% | -42.2 | 0% | 0 | 0 |

## Branch sizes on the omitted-loop backbone (graph 15), median fitted vs truth

Truth columns are the true effective trajectory: size at the young end, old end, its harmonic mean
(accumulated drift ⇔ duration / harmonic) and arithmetic mean. For b the truth is the composite of the
leaf, the loop (1 / Σ p_k²/N_k) and anc_b. Fitted = median of posterior means across subsamples, at the
best-ELBO event order of g15 on each subsample (the same fits the figures use).

### growth

| branch | young | old | harmonic | arithmetic | true shared | true sep IBD | true sep SNP | hapibd shared | hapibd sep IBD | hapibd sep SNP |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| a | 150,000 | 15,000 | 38,376 | 58,630 | 39,809 | 40,121 | 17,192 | 37,730 | 39,473 | 16,632 |
| b | 60,000 | 15,000 | 20,776 | 23,320 | 17,595 | 19,234 | 18,916 | 18,455 | 19,089 | 18,488 |
| c | 15,000 | 15,000 | 15,000 | 15,000 | 15,861 | 14,306 | 17,088 | 15,961 | 14,514 | 17,511 |
| b.1 | 15,000 | 15,000 | 15,000 | 15,000 | 14,428 | 15,960 | 19,340 | 15,522 | 16,487 | 18,168 |
| b.2 | 5,000 | 5,000 | 5,000 | 5,000 | 7,559 | 25,750 | 16,915 | 7,840 | 29,663 | 17,488 |
| left | 15,000 | 15,000 | 15,000 | 15,000 | 12,436 | 10,746 | 18,160 | 13,739 | 12,297 | 17,232 |
| right | 15,000 | 15,000 | 15,000 | 15,000 | 9,691 | 24,675 | 16,542 | 9,473 | 30,077 | 17,717 |
| root | 15,000 | 15,000 | 15,000 | 15,000 | 13,639 | 15,359 | 16,795 | 13,143 | 16,788 | 17,107 |

### constant

| branch | young | old | harmonic | arithmetic | true shared | true sep IBD | true sep SNP | hapibd shared | hapibd sep IBD | hapibd sep SNP |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| a | 15,000 | 15,000 | 15,000 | 15,000 | 14,327 | 14,108 | 15,054 | 14,403 | 14,453 | 15,268 |
| b | 15,000 | 15,000 | 18,462 | 19,500 | 15,075 | 16,614 | 15,537 | 15,724 | 16,219 | 15,954 |
| c | 15,000 | 15,000 | 15,000 | 15,000 | 15,676 | 14,981 | 16,800 | 16,104 | 16,067 | 14,290 |
| b.1 | 15,000 | 15,000 | 15,000 | 15,000 | 14,698 | 15,457 | 15,760 | 15,997 | 16,097 | 15,881 |
| b.2 | 15,000 | 15,000 | 15,000 | 15,000 | 16,132 | 19,852 | 17,256 | 15,445 | 16,811 | 15,420 |
| left | 15,000 | 15,000 | 15,000 | 15,000 | 14,550 | 13,930 | 15,535 | 16,147 | 15,672 | 15,421 |
| right | 15,000 | 15,000 | 15,000 | 15,000 | 15,307 | 19,552 | 17,029 | 15,598 | 16,309 | 14,367 |
| root | 15,000 | 15,000 | 15,000 | 15,000 | 15,236 | 15,683 | 15,441 | 15,529 | 16,258 | 14,871 |

## Residual tilt on the backbone (Pearson residual per cM), median across subsamples

| Scenario | IBD | Model | a–a | a–b | a–c | b–b | b–c | c–c |
|---|---|---|---:|---:|---:|---:|---:|---:|
| growth | true | poisson_shared | -0.018 | -0.001 | +0.007 | +0.005 | -0.006 | -0.019 |
| growth | true | poisson_separate | -0.026 | -0.022 | +0.001 | -0.014 | +0.043 | -0.031 |
| growth | hapibd | poisson_shared | -0.013 | -0.000 | +0.007 | -0.007 | -0.004 | -0.018 |
| growth | hapibd | poisson_separate | -0.033 | -0.020 | +0.002 | -0.012 | +0.050 | -0.036 |
| constant | true | poisson_shared | -0.000 | +0.004 | +0.003 | +0.002 | +0.004 | +0.022 |
| constant | true | poisson_separate | +0.012 | -0.019 | +0.002 | +0.008 | +0.024 | +0.012 |
| constant | hapibd | poisson_shared | -0.007 | +0.000 | +0.002 | +0.007 | +0.004 | +0.022 |
| constant | hapibd | poisson_separate | +0.006 | -0.004 | +0.002 | +0.004 | +0.003 | +0.016 |

Figures: `figures/` — topology_recovery, evidence_benefit, spectrum_fit_*, ne_trajectory_*, ne_ibd_vs_snp_*,
residual_tilt_*, component_likelihood_*, parameters_backbone_*, topology_generating_*.
