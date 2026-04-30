data {
    // Topology block
    int<lower=0> n_leaves;
    int<lower=0> n_events;
    int<lower=0> n_nodes;
    int<lower=0> n_admixture;

    array[n_events] matrix[n_nodes, n_nodes] migration_matrices;
    array[n_admixture, 4] int<lower=0> admixture_map;

    array[n_admixture] int<lower=1, upper=n_events+1> admixture_indices;
    array[n_events - n_admixture] int<lower=1, upper=n_events> fixed_indices;
    array[n_events - n_admixture] int<lower=1, upper=n_events+1> fixed_indices_shifted;

    // Leaf pairs
    int<lower=0> n_leaf_pairs;
    array[n_leaf_pairs] int<lower=1, upper=n_leaves> pair_i;
    array[n_leaf_pairs] int<lower=1, upper=n_leaves> pair_j;

    // Observation block
    matrix[n_leaves, n_leaves] w_hat;
    matrix<lower=0>[n_leaves, n_leaves] w_se;
}

parameters {
    vector<lower=1>[n_events] times;
    array[n_admixture] real<lower=0, upper=1> admixture_fractions;
    vector<lower=1000>[n_nodes] Ne;
}

transformed parameters {
    vector[n_events] cumulative_times = cumulative_sum(times);
    matrix[n_nodes, n_nodes] I = diag_matrix(rep_vector(1.0, n_nodes));

    // ================================================================
    // 1. Build parameter_migration_matrices
    // ================================================================
    array[n_events + 1] matrix[n_nodes, n_nodes] parameter_migration_matrices;
    parameter_migration_matrices[1] = I;

    for (k in 1:(n_events - n_admixture)) {
        parameter_migration_matrices[fixed_indices_shifted[k]] =
            migration_matrices[fixed_indices[k]];
    }

    for (i in 1:n_admixture) {
        int src  = admixture_map[i, 2]+1;
        int tgt1 = admixture_map[i, 3]+1;
        int tgt2 = admixture_map[i, 4]+1;

        matrix[n_nodes, n_nodes] admix_mat = I;
        admix_mat[src, tgt1] =  admixture_fractions[i];
        admix_mat[src, tgt2] =  1.0 - admixture_fractions[i];
        admix_mat[src,  src] =  0.0;

        parameter_migration_matrices[admixture_indices[i]] = admix_mat;
    }

    // ================================================================
    // 2. SNP: Brownian drift accumulation with per-node Ne
    // ================================================================
    matrix[n_leaves, n_leaves] snp_mean = rep_matrix(0.0, n_leaves, n_leaves);
    for (p in 1:n_leaf_pairs) {
        int li = pair_i[p];
        int lj = pair_j[p];

        matrix[n_nodes, n_nodes] W = rep_matrix(0.0, n_nodes, n_nodes);
        W[li, lj] = 1.0;

        for (e in 1:n_events) {
            if (e > 1) {
                matrix[n_nodes, n_nodes] M = parameter_migration_matrices[e];
                W = M' * W * M;
            }

            real duration;
            if (e == 1) {
                duration = cumulative_times[1];
            } else {
                duration = cumulative_times[e] - cumulative_times[e - 1];
            }

            for (a in 1:n_nodes) {
                real w_diag = W[a, a];
                if (w_diag > 1e-20) {
                    snp_mean[li, lj] += w_diag * duration / Ne[a];
                }
            }
        }

        if (li != lj) {
            snp_mean[lj, li] = snp_mean[li, lj];
        }
    }

    // ================================================================
    // 3. Double-centering
    // ================================================================
    matrix[n_leaves, n_leaves] W_centered;
    {
        vector[n_leaves] row_means;
        real grand_mean;

        for (i in 1:n_leaves) {
            row_means[i] = mean(snp_mean[i]);
        }
        grand_mean = mean(row_means);

        W_centered = snp_mean - rep_matrix(row_means, n_leaves)
              - rep_matrix(row_means', n_leaves)
              + grand_mean;
    }
}

model {
    times ~ exponential(1.0 / 100);
    admixture_fractions ~ uniform(0, 1);
    Ne ~ gamma(2, 0.0002);

    for (i in 1:n_leaves) {
        for (j in i:n_leaves) {
            w_hat[i, j] ~ normal(W_centered[i, j], w_se[i, j]);
        }
    }
}
