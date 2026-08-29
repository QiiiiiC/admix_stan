functions {
    real int_p_L(real N, real t1, real t2, real u, real v) {
        if (t2 - t1 < 1e-10) return 0.0;  // zero-length epoch

        real cu = -(1.0/N + u/50.0);
        real cv = -(1.0/N + v/50.0);
        
        // Push exp(t1/N) inside the calculation to avoid Inf * 0 = NaN
        real eu1 = exp(-u * t1 / 50.0);
        real eu2 = exp((t1 - t2) / N - u * t2 / 50.0);
        real ev1 = exp(-v * t1 / 50.0);
        real ev2 = exp((t1 - t2) / N - v * t2 / 50.0);
        
        real cu2 = cu * cu;
        real cv2 = cv * cv;

        real k1 = u/(50.0*N) * (eu2 * (cu*t2 - 1.0)/cu2 - eu1 * (cu*t1 - 1.0)/cu2);
        real k2 = v/(50.0*N) * (ev2 * (cv*t2 - 1.0)/cv2 - ev1 * (cv*t1 - 1.0)/cv2);
        
        real k3 = 1.0/N * (eu2/cu - eu1/cu);
        real k4 = 1.0/N * (ev2/cv - ev1/cv);

        return k1 - k2 + k3 - k4;
    }

    real int_p_LL(real N, real t1, real t2, real u, real v){
        if (t2 - t1 < 1e-10) return 0.0;  // zero-length epoch

        real cu = -(1.0/N + u/50.0);
        real cv = -(1.0/N + v/50.0);
        
        real eu1 = exp(-u * t1 / 50.0);
        real eu2 = exp((t1 - t2) / N - u * t2 / 50.0);
        real ev1 = exp(-v * t1 / 50.0);
        real ev2 = exp((t1 - t2) / N - v * t2 / 50.0); 

        
        real k1 = 1.0/(50.0*N) * (eu2*(cu*t2 - 1.0)/pow(cu, 2) - eu1*(cu*t1 - 1.0)/pow(cu, 2));
        real k2 = 1.0/(50.0*N) * (ev2*(cv*t2 - 1.0)/pow(cv, 2) - ev1*(cv*t1 - 1.0)/pow(cv, 2));

        real result = k1 - k2;
        return result;
    }
}

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

    // Fixed parameter block
    int<lower=0> n_bins;
    array[n_bins, 2] real<lower=0> bin_length;
    real<lower=0> T_max;
    real<lower=0> cm;

    // Precomputed leaf pair indices (1-based, pair_i[p] <= pair_j[p])
    int<lower=0> n_leaf_pairs;
    array[n_leaf_pairs] int<lower=1, upper=n_leaves> pair_i;
    array[n_leaf_pairs] int<lower=1, upper=n_leaves> pair_j;

    // Number of haploid samples per population (leaf), used for pair counts
    array[n_leaves] int<lower=0> n_samples;

    // Observation block
    array[n_bins] matrix<lower=0>[n_leaves, n_leaves] ibd_hat;
    array[n_bins] matrix<lower=0>[n_leaves, n_leaves] ibd_se;
}

parameters {
    vector<lower=1>[n_events] times;
    array[n_admixture] real<lower=0, upper=1> admixture_fractions;
    vector<lower=0>[n_events + 1] Ne_epoch;  // one shared Ne per epoch
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
    // 2. Compute IBD fractions via joint path tracking
    //
    // W[a, b] = joint probability weight for (lineage_i at a, lineage_j at b)
    //           accounting for no-coalescence so far.
    //
    // The key question: what does W carry exactly?
    //
    // We define W so that the IBD contribution from epoch e is:
    //   sum_a W[a,a] * int_p_L(N_a, t_start, t_end, u, v)
    //
    // Since int_p_L(N, t1, t2, u, v) is the FULL conditional IBD
    // integral (given lineages are co-located in pop N at t1), the
    // weight W[a,a] must carry:
    //   - path probability (admixture fractions)
    //   - survival (no coalescence) in all PREVIOUS epochs
    //
    // Survival in epoch with duration d in population of size N is
    // S = exp(-d/N). But int_p_L already "knows" the epoch starts at
    // t_start — it integrates (1/N)exp(-(t - t_start)/N).
    //
    // So W[a,a] after epoch e must be:
    //   (path weight) * prod_{prev epochs} S(N_k, duration_k)
    //
    // where S(N, d) = exp(-d/N) = probability of no coalescence in
    // a population of size N for d generations.
    //
    // After accumulating IBD from epoch e, we multiply:
    //   W[a,a] *= S(N_a, duration_e) = exp(-duration_e / N_a)
    //
    // For off-diagonal W[a,b] (a != b): lineages can't coalesce when
    // in different nodes, so no survival factor is needed. W[a,b]
    // passes through unchanged.
    // ================================================================

    array[n_bins] matrix[n_leaves, n_leaves] ibd_fraction;
    array[n_bins] matrix[n_leaves, n_leaves] ibd_number;
    for (b in 1:n_bins) {
        ibd_fraction[b] = rep_matrix(0.0, n_leaves, n_leaves);
        ibd_number[b] = rep_matrix(0.0, n_leaves, n_leaves);
    }

    for (p in 1:n_leaf_pairs) {
        int li = pair_i[p];
        int lj = pair_j[p];

        // Initialize: both lineages at their starting leaf nodes
        matrix[n_nodes, n_nodes] W = rep_matrix(0.0, n_nodes, n_nodes);
        W[li, lj] = 1.0;

        for (e in 1:(n_events + 1)) {
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
            } else if (e == n_events + 1) {
                t_start = cumulative_times[n_events];
                t_end = T_max;
            } else {
                t_start = cumulative_times[e - 1];
                t_end = cumulative_times[e];
            }

            real duration = t_end - t_start;

            // --- IBD accumulation + survival update ---
            for (a in 1:n_nodes) {
                real w_diag = W[a, a];

                if (w_diag > 1e-20) {
                    // IBD: W[a,a] carries path weights * prior survival.
                    // int_p_L gives the conditional IBD integral for this epoch.
                    for (b in 1:n_bins) {
                        ibd_fraction[b][li, lj] += w_diag *
                            int_p_L(Ne_epoch[e], t_start, t_end,
                                    bin_length[b, 1], bin_length[b, 2]);
                        ibd_number[b][li, lj] += w_diag * int_p_LL(Ne_epoch[e], t_start, t_end,
                                    bin_length[b, 1], bin_length[b, 2]);
                    }

                    // Survival: no coalescence in this epoch
                    W[a, a] = w_diag * exp(-duration / Ne_epoch[e]);
                }
            }
        }

        // Fill symmetric entry
        if (li != lj) {
            for (b in 1:n_bins) {
                ibd_fraction[b][lj, li] = ibd_fraction[b][li, lj];
                ibd_number[b][lj, li] = ibd_number[b][li, lj];
            }
        }
    }
}

model {
    times ~ exponential(0.01);
    admixture_fractions ~ beta(1, 1);
    // Gamma prior: mean 15000, sd 7500 (shape 4 => coefficient of variation 0.5)
    Ne_epoch ~ gamma(4.0, 4.0 / 15000.0);

    for (i in 1:n_leaves) {
        for (j in i:n_leaves) {
            for (b in 1:n_bins) {

                    // MASKED: zero-observation Poisson penalty removed; only the
                    // nonzero (observed) IBD cells enter the likelihood.
                    if(ibd_hat[b][i, j] >0){
                        target += normal_lpdf(ibd_hat[b][i,j] |ibd_fraction[b][i, j], ibd_se[b][i, j]);
                    }
                
            }
        }
    }
}

