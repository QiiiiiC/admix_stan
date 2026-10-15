// Relaxed graph, E-step: theta is DATA; the random object is the slot logits y.
//   y ~ Gumbel(0, 1) iid          =>  argmax_c y_kc ~ Uniform(C_k)      (Gumbel-max, eq. 15)
//   w_k = softmax(y_k / lambda)   soft slot weights; M_k = sum_c w_kc M_k^(c)   (eq. 16)
//   target += l(theta, M(w))
// As lambda -> 0 the posterior of argmax y tends to p(v | x, theta)   (eq. 17-18).
functions {
#include graph.stanfunctions
}
data {
    int<lower=1> L;
    int<lower=1> n_merge;
    int<lower=1> n_pairs;
    array[n_pairs, 2] int<lower=1> pairs;
    real<lower=0> lambda;
    vector<lower=1>[n_merge + 2] times;
    real<lower=0, upper=1> f;
    real<lower=1> N;
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
transformed data {
    vector[n_merge + 2] t = cumulative_sum(times);
}
parameters {
    vector[L] y_admix;
    array[n_merge] vector[n_pairs] y_merge;
}
transformed parameters {
    vector[L] w_admix = softmax(y_admix / lambda);
    array[n_merge] vector[n_pairs] w_merge;
    for (k in 1:n_merge) w_merge[k] = softmax(y_merge[k] / lambda);
    real loglik = sum(graph_loglik(slot_matrices(w_admix, w_merge, f, pairs, L), t, N, bin_length, T_max, cm,
                                   n_samples, ibd_count, w_hat, w_se, use_snp));
}
model {
    y_admix ~ gumbel(0, 1);
    for (k in 1:n_merge) y_merge[k] ~ gumbel(0, 1);
    target += loglik;
}
generated quantities {
    int hard_admix = sort_indices_desc(y_admix)[1];
    array[n_merge] int hard_merge;
    for (k in 1:n_merge) hard_merge[k] = sort_indices_desc(y_merge[k])[1];
}
