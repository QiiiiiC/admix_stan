// Discrete graphs, theta = (event-time increments, f, N) fixed unknown constants.
// target = log p(x | theta) = log sum_v p(v) exp l(theta, v)            (MATHS.md eq. 12)
// Maximising it (optimize, no Jacobian) is the EM / type-II ML estimate of theta; with
// one vertex it is that graph's own maximum likelihood (the profile score).
// Generated quantities: p(v | x, theta) (eq. 13) and every vertex's log-likelihood.
functions {
#include graph.stanfunctions
}
data {
    int<lower=1> L;
    int<lower=1> n_merge;                          // merge slots between admixture and root
    int<lower=1> n_pairs;
    array[n_pairs, 2] int<lower=1> pairs;
    int<lower=1> V;                                 // vertices (discrete graphs) summed over
    array[V] int<lower=1, upper=L> vertex_admix;
    array[V, n_merge] int<lower=1, upper=n_pairs> vertex_merge;
    vector[V] log_prior;
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
    int K = n_merge + 2;
}
parameters {
    vector<lower=1>[K] times;
    real<lower=0, upper=1> f;
    real<lower=1> N;
}
transformed parameters {
    vector[K] t = cumulative_sum(times);
    vector[V] loglik;
    for (v in 1:V) {
        vector[L] wa = one_hot_vector(L, vertex_admix[v]);
        array[n_merge] vector[n_pairs] wm;
        for (k in 1:n_merge) wm[k] = one_hot_vector(n_pairs, vertex_merge[v, k]);
        loglik[v] = sum(graph_loglik(slot_matrices(wa, wm, f, pairs, L), t, N, bin_length, T_max, cm,
                                     n_samples, ibd_count, w_hat, w_se, use_snp));
    }
}
model {
    target += log_sum_exp(log_prior + loglik);
}
generated quantities {
    vector[V] posterior = softmax(log_prior + loglik);
}
