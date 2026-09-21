# Evidence-estimator check on the winning backbone (g15_o1, growth, true IBD)

Two quick alternatives to the Pathfinder importance-sampling logZ, on the same 10 subsamples and the same two models. Laplace: best of four L-BFGS starts (three dispersed, one at Pathfinder's best draw), CmdStan `laplace` draws (4,000), logZ = lp(mode) + d/2 log 2π + ½ log|Σ| with Σ the draw covariance in unconstrained space; the same Gaussian used as an IS proposal. Student-t IS: ν = 4, scale 1.5× the covariance, 4,000 draws, centred on the pooled Pathfinder draws or on the Laplace mode; log p evaluated with CmdStan `log_prob` (jacobian on). Draws where the model throws (overflow far in the tail) are dropped and counted.

## Pareto k and ESS per estimator (median [min, max] across the 10 subsamples)

| model | estimator | Pareto k | ESS of 4,000 (12,000 for Pathfinder) | dropped draws |
|---|---|---|---|---|
| shared | Pathfinder IS | 1.62 [1.52, 1.83] | 3 [1, 14] | 0 [0, 0] |
| shared | Laplace-Gaussian IS | 1.16 [0.94, 7.86] | 8 [1, 30] | 0 [0, 0] |
| shared | Student-t @Pathfinder | 2.16 [1.82, 5.56] | 3 [1, 5] | 2 [0, 59] |
| shared | Student-t @Laplace | 5.45 [4.60, 11.83] | 4 [1, 12] | 36 [22, 61] |
| separate | Pathfinder IS | 2.14 [1.73, 2.37] | 3 [2, 9] | 0 [0, 0] |
| separate | Laplace-Gaussian IS | 5.46 [4.80, 6.15] | 9 [1, 17] | 0 [0, 0] |
| separate | Student-t @Pathfinder | 32.87 [4.11, 44.70] | 1 [1, 2] | 370 [3, 4000] |
| separate | Student-t @Laplace | 16.49 [12.07, 20.87] | 3 [1, 6] | 15 [7, 20] |

None of the estimators reaches k ≤ 0.7 on any subsample. Heavier tails make k worse, not better.

## logZ per subsample (nats)

| model | rep | Pathfinder ELBO | Pathfinder IS logZ | Laplace (analytic) | Laplace-Gaussian IS | Student-t IS @Pathfinder | Student-t IS @Laplace mode |
|---|---|---|---|---|---|---|---|
| shared | 0 | -206.1 | -191.7 | -179.7 | -181.4 | -186.7 | -182.2 |
| shared | 1 | -209.9 | -192.9 | -180.7 | -182.9 | -187.6 | -182.0 |
| shared | 2 | -213.1 | -192.4 | -176.8 | -179.8 | -181.7 | -179.0 |
| shared | 3 | -216.1 | -197.7 | -180.7 | -182.8 | -187.6 | -182.4 |
| shared | 4 | -212.0 | -191.6 | -175.7 | -177.7 | -183.3 | -178.1 |
| shared | 5 | -208.9 | -194.4 | -179.1 | -180.6 | -185.8 | -181.6 |
| shared | 6 | -215.2 | -193.3 | -182.2 | 75.3 | -187.1 | -185.7 |
| shared | 7 | -204.9 | -187.4 | -176.2 | -178.4 | -181.8 | -176.2 |
| shared | 8 | -231.0 | -206.3 | -191.5 | -193.0 | -199.9 | -193.5 |
| shared | 9 | -221.7 | -205.8 | -191.6 | -193.9 | -196.0 | -193.1 |
| separate | 0 | -278.6 | -203.6 | -180.4 | -186.0 | -210.8 | -186.0 |
| separate | 1 | -404.6 | -206.4 | -180.3 | -185.3 | n/a | -185.5 |
| separate | 2 | -391.6 | -205.5 | -177.6 | -183.0 | n/a | -183.4 |
| separate | 3 | -413.8 | -212.6 | -182.0 | -187.2 | n/a | -187.0 |
| separate | 4 | -279.8 | -209.7 | -176.5 | -181.2 | -217.7 | -176.7 |
| separate | 5 | -287.4 | -209.7 | -178.2 | -183.8 | -207.9 | -184.1 |
| separate | 6 | -243.9 | -210.5 | -182.5 | -186.2 | -198.3 | -188.1 |
| separate | 7 | -271.3 | -202.7 | -178.8 | -183.3 | -197.5 | -183.1 |
| separate | 8 | -292.9 | -227.2 | -191.5 | -198.0 | -224.8 | -196.4 |
| separate | 9 | -288.1 | -228.0 | -192.5 | -197.4 | -220.1 | -198.6 |

## Paired separate − shared (nats; > 0 favours separate Ne)

| estimator | mean | subsamples > 0 |
|---|---|---|
| Pathfinder ELBO | -101.3 | 0/10 |
| Pathfinder IS logZ | -16.2 | 0/10 |
| Laplace (analytic) | -0.6 | 2/10 |
| Laplace-Gaussian IS | -29.6 | 0/10 |
| Student-t IS @Pathfinder | -22.4 | 0/7 |
| Student-t IS @Laplace mode | -3.5 | 1/10 |

## Why the Gaussian and Student-t proposals fail

Eigen-decomposing the Laplace covariance in unconstrained space (shared model) shows the problem is the *shape* of the posterior, not its tails. Laplace draws within one standard deviation of the mode along the widest axis have a mean log-posterior 478 nats (subsample 0) and 4,056 nats (subsample 6) below the mode; a quadratic log-posterior would give about −8. The Hessian at the mode says that direction is flat (sd ≈ 1), yet Pathfinder's draws along it have sd 0.09 (subsample 0): the posterior is a narrow ridge that curves away from the mode's tangent. The widest axes load on the deep-node sizes (Ne_raw of root, n1, n2 with log sigma) and on the shallow (b admixture time, Ne of b) pair. A Gaussian centred at the mode puts most of its mass in empty space, so the analytic Laplace logZ is biased upward (it is 12–30 nats above every sampling-based estimate), and a heavier-tailed proposal only sends more draws into empty space while the rare draw that lands on the ridge takes all the weight.

Consequence: any estimator built on a single ellipsoid (Pathfinder, Laplace, Student-t IS) cannot deliver a converged logZ for this model. The ELBO remains the only stable Pathfinder statistic. Converged evidence needs draws that follow the ridge: NUTS with bridge sampling, sequential Monte Carlo, or a reparameterisation in drift coordinates that straightens it.

## What this does to the shared-vs-separate conclusion

The Pathfinder ELBO says separate Ne loses by about 100 nats on this graph; the Laplace evidence says the two models are tied (mean −0.6, two subsamples positive). The difference is the approximation gap, not the data. On the separate model the Pathfinder draws sit on average 100–240 nats below the posterior mode (ELBO −244 to −414 against a mode log-density near −175), whereas on the shared model the gap is about 35 nats. Pathfinder's Gaussian mixture covers the 27-dimensional separate posterior far worse than the 16-dimensional shared one, and the ELBO charges that to the model. The separate model's mode log-density is 7 nats below the shared one's, which is the cost of eleven extra prior densities evaluated at their centres, and the Laplace volume term gives that back.

So "separate − shared < 0 in 80/80 comparisons" is a statement about the ELBO, and the ELBO is not a fair judge between models of different dimension when its gap differs by model. Neither estimator here is converged, so the honest summary is: the current evidence does not favour separate Ne, and also does not show it losing. Settling it needs a ridge-following estimator (NUTS plus bridge sampling, or SMC) or the drift reparameterisation.
