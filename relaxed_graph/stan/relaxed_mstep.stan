// Relaxed graph, M-step: theta maximises the Monte Carlo E-step objective
//   Q(theta) = (1/S) sum_s l(theta, M(w_s))                               (eq. 19)
// over S soft-weight draws w_s from the E-step posterior.  Run with optimize (no Jacobian).
functions {
#include graph.stanfunctions
}
data {
    int<lower=1> L;
    int<lower=1> n_merge;
    int<lower=1> n_pairs;
    array[n_pairs, 2] int<lower=1> pairs;
    int<lower=1> S;
    array[S] vector[L] w_admix;
    array[S, n_merge] vector[n_pairs] w_merge;
    int<lower=1> n_bins;
    array[n_bins, 2] real<lower=0> bin_length;
    real<lower=0> T_max;
    real<lower=0> cm;
    array[L] int<lower=0> n_samples;
    array[n_bins, L, L] int<lower=0> ibd_count;
    matrix[L, L] w_hat;
    matrix<lower=0>[L, L] w_se;
    int<lower=0, upper=1> use_snp;
}
parameters {
    vector<lower=1>[n_merge + 2] times;
    real<lower=0, upper=1> f;
    real<lower=1> N;
}
model {
    vector[n_merge + 2] t = cumulative_sum(times);
    real q = 0;
    for (s in 1:S)
        q += sum(graph_loglik(slot_matrices(w_admix[s], w_merge[s], f, pairs, L), t, N, bin_length, T_max, cm,
                              n_samples, ibd_count, w_hat, w_se, use_snp));
    target += q / S;
}
