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

    // Nsmooth: tree structure for the log-Ne random walk (parents have a
    // higher node index than their children).
    array[n_nodes] int<lower=0> ne_parent;       // parent node (1-based); 0 = root/admix
    array[n_nodes] int<lower=0> ne_admix_idx;    // admixture index if admix-src node; 0 else
    array[n_nodes] int<lower=0> ne_start_event;  // creating event (1-based); 0 = leaf (t=0)

    // Fixed parameter block
    int<lower=0> n_leaf_pairs;
    array[n_leaf_pairs] int<lower=1, upper=n_leaves> pair_i;
    array[n_leaf_pairs] int<lower=1, upper=n_leaves> pair_j;

    // Observation block
    matrix[n_leaves, n_leaves] w_hat;
    matrix<lower=0>[n_leaves, n_leaves] w_se;
}

parameters{
    vector<lower = 1>[n_events] times;
    array[n_admixture] real<lower = 0, upper = 1> admixture_fractions;

    // Hierarchical log-normal for Ne, non-centered, as a tree-structured
    // Gaussian random walk -- identical construction to ibd_model_Nsmooth so the
    // two likelihoods are compared under the SAME Ne parametrisation.
    real mu_log;
    real<lower=0> sigma_log;
    real<lower=0> tau;            // log-Ne random-walk step scale (per t_ref gen)
    vector[n_nodes] Ne_raw;       // std-normal RW increments (non-centered)
}

transformed parameters{
    vector[n_events] cumulative_times = cumulative_sum(times);

    // === Nsmooth: tree-structured Gaussian random walk on log-Ne ===
    // Brownian motion in TIME along the tree: each branch is centred on the one
    // it flows into going backwards (its parent), with variance tau^2*dt/t_ref
    // for a branch of duration dt.  sqrt(dt) makes the prior invariant to how
    // finely the graph is subdivided.  Filled root->leaves by descending index.
    real t_ref = 100.0;
    vector[n_nodes] node_t;
    for (a in 1:n_nodes)
        node_t[a] = ne_start_event[a] == 0 ? 0.0 : cumulative_times[ne_start_event[a]];

    vector[n_nodes] log_Ne;
    for (ii in 1:n_nodes) {
        int a = n_nodes - ii + 1;
        if (ne_parent[a] == 0 && ne_admix_idx[a] == 0) {
            log_Ne[a] = mu_log + sigma_log * Ne_raw[a];            // root anchor
        } else if (ne_admix_idx[a] > 0) {
            int i  = ne_admix_idx[a];
            int p1 = admixture_map[i, 3] + 1;
            int p2 = admixture_map[i, 4] + 1;
            real fr = admixture_fractions[i];
            real pmean = fr * log_Ne[p1] + (1.0 - fr) * log_Ne[p2];
            real dt = node_t[p1] - node_t[a];
            log_Ne[a] = pmean + tau * sqrt(dt / t_ref + 1e-9) * Ne_raw[a];
        } else {
            int p  = ne_parent[a];
            real dt = node_t[p] - node_t[a];
            log_Ne[a] = log_Ne[p] + tau * sqrt(dt / t_ref + 1e-9) * Ne_raw[a];
        }
    }
    vector<lower=0>[n_nodes] Ne = exp(log_Ne);

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
    
    matrix[n_leaves, n_leaves] snp_mean = rep_matrix(0.0, n_leaves, n_leaves);
    for (p in 1:n_leaf_pairs) {
        int li = pair_i[p];
        int lj = pair_j[p];

        // Initialize: both lineages at their starting leaf nodes
        matrix[n_nodes, n_nodes] W = rep_matrix(0.0, n_nodes, n_nodes);
        W[li, lj] = 1.0;

        for (e in 1:(n_events)) {
            // --- Migration event before epoch e ---
            if (e > 1) {
                matrix[n_nodes, n_nodes] M = parameter_migration_matrices[e];
                W = M' * W * M;
            }

            // --- Epoch boundaries ---
            real t_start;
            real t_end;
            if (e == 1) {
                t_start = 0.0;
                t_end = cumulative_times[1];
            } else {
                t_start = cumulative_times[e - 1];
                t_end = cumulative_times[e];
            }

            real duration = t_end - t_start;

            // --- SNP drift accumulation (Brownian, no coalescence tracking) ---
            for (a in 1:n_nodes) {
                real w_diag = W[a, a];

                if (w_diag > 1e-20) {
                    snp_mean[li, lj] += w_diag * duration / Ne[a];
                }
            }
        }

        // Fill symmetric entry
        if (li != lj) {
            snp_mean[lj, li] = snp_mean[li, lj];
        }
        
    }


    matrix[n_leaves, n_leaves] W_centered;
    {
        vector[n_leaves] row_means;
        real grand_mean;

        // Calculate Row Means
        for (i in 1:n_leaves) {
            row_means[i] = mean(snp_mean[i]);
        }

        // Calculate Grand Mean
        grand_mean = mean(row_means);

        W_centered = snp_mean - rep_matrix(row_means, n_leaves)
              - rep_matrix(row_means', n_leaves)
              + grand_mean;
    }
}

model {
    times ~ exponential(1.0/100);
    admixture_fractions ~ beta(1.0,1.0);
    mu_log    ~ normal(log(15000), 0.53);  // log-scale location; median Ne = 15000
    // 0.53 = sqrt(trigamma(4)), the log-scale sd of the gamma(4, 4/15000) prior
    // that Nfixed puts on effective_N -- so Nsmooth and Nfixed now express the
    // SAME prior belief about the LEVEL of Ne and differ only in whether it may
    // vary across the graph.  At the old 0.25 the two were not comparable: +-28%
    // at 1 sigma made Ne = 3,300 a ~6-sigma excursion, so SNP-only-Nsmooth on the
    // AFR-EUR two-leaf test reported the prior (Ne = 10,100, g ratio 0.35) while
    // SNP-only-Nfixed reached the data (Ne = 3,323, g ratio 1.00).
    sigma_log ~ normal(0, 0.3);            // half-normal; per-node spread
    tau       ~ normal(0, 0.3);            // half-normal; log-Ne RW step scale
    Ne_raw    ~ std_normal();

    for (i in 1:n_leaves) {
        for (j in i:n_leaves) {
            w_hat[i, j] ~ normal(W_centered[i, j], w_se[i, j]);
        }
    }
}

generated quantities {
    real lp_snp;
    real chi2_snp;
    int n_snp_obs;

    lp_snp = 0;
    chi2_snp = 0;
    for (i in 1:n_leaves) {
        for (j in i:n_leaves) {
            lp_snp += normal_lpdf(w_hat[i,j] | W_centered[i,j], w_se[i,j]);
            chi2_snp += square((w_hat[i,j] - W_centered[i,j]) / w_se[i,j]);
        }
    }
    n_snp_obs = n_leaves * (n_leaves + 1) %/% 2;
}
