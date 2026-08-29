functions {
    real int_p_L(real N, real t1, real t2, real u, real v) {
        if (t2 - t1 < 1e-10) return 0.0;

        real cu = -(1.0/N + u/50.0);
        real cv = -(1.0/N + v/50.0);

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
        if (t2 - t1 < 1e-10) return 0.0;

        real cu = -(1.0/N + u/50.0);
        real cv = -(1.0/N + v/50.0);

        real eu1 = exp(-u * t1 / 50.0);
        real eu2 = exp((t1 - t2) / N - u * t2 / 50.0);
        real ev1 = exp(-v * t1 / 50.0);
        real ev2 = exp((t1 - t2) / N - v * t2 / 50.0);

        real k1 = 1.0/(50.0*N) * (eu2*(cu*t2 - 1.0)/pow(cu, 2) - eu1*(cu*t1 - 1.0)/pow(cu, 2));
        real k2 = 1.0/(50.0*N) * (ev2*(cv*t2 - 1.0)/pow(cv, 2) - ev1*(cv*t1 - 1.0)/pow(cv, 2));

        return k1 - k2;
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

    // Nsmooth: tree structure for the log-Ne random walk (parents have a
    // higher node index than their children).
    array[n_nodes] int<lower=0> ne_parent;       // parent node (1-based); 0 = root/admix
    array[n_nodes] int<lower=0> ne_admix_idx;    // admixture index if admix-src node; 0 else
    array[n_nodes] int<lower=0> ne_start_event;  // creating event (1-based); 0 = leaf (t=0)

    // Fixed parameter block
    int<lower=0> n_bins;
    array[n_bins, 2] real<lower=0> bin_length;
    real<lower=0> T_max;
    real<lower=0> cm;

    // Precomputed leaf pair indices
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

    // Hierarchical log-normal prior for effective population sizes (non-centered):
    //   Ne[a] = exp(mu_log + sigma_log * Ne_raw[a]),  Ne_raw[a] ~ N(0, 1).
    // mu_log = mean-of-log (normal prior), sigma_log = sd-of-log (half-normal).
    // Ne is built in the transformed parameters block below.
    real mu_log;
    real<lower=0> sigma_log;
    real<lower=0> tau;            // log-Ne random-walk step scale (per t_ref gen)
    vector[n_nodes] Ne_raw;       // std-normal RW increments (non-centered)

}

transformed parameters {

    vector[n_events] cumulative_times = cumulative_sum(times);

    // === Nsmooth: tree-structured Gaussian random walk on log-Ne ===
    // Each branch's Ne is centred on the branch it flows INTO going backwards in
    // time (its parent); admixture branches are centred on the fraction-weighted
    // mean of their two source branches.  Variance = tau^2 * dt / t_ref, dt =
    // branch duration.  Parents have higher node index than children, so we fill
    // log_Ne from the root (highest index) down to the leaves.
    real t_ref = 100.0;                       // reference branch length (generations)
    vector[n_nodes] node_t;                   // branch start time (younger end)
    for (a in 1:n_nodes)
        node_t[a] = ne_start_event[a] == 0 ? 0.0 : cumulative_times[ne_start_event[a]];

    vector[n_nodes] log_Ne;
    for (ii in 1:n_nodes) {
        int a = n_nodes - ii + 1;             // descending index: parents before children
        if (ne_parent[a] == 0 && ne_admix_idx[a] == 0) {
            log_Ne[a] = mu_log + sigma_log * Ne_raw[a];            // root anchor
        } else if (ne_admix_idx[a] > 0) {
            int i  = ne_admix_idx[a];
            int p1 = admixture_map[i, 3] + 1;                      // source with fraction f
            int p2 = admixture_map[i, 4] + 1;                      // source with fraction 1-f
            real f = admixture_fractions[i];
            real pmean = f * log_Ne[p1] + (1.0 - f) * log_Ne[p2];
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

    // ================================================================
    // 2. Compute IBD fractions via joint path tracking
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

        matrix[n_nodes, n_nodes] W = rep_matrix(0.0, n_nodes, n_nodes);
        W[li, lj] = 1.0;

        for (e in 1:(n_events + 1)) {
            if (e > 1) {
                matrix[n_nodes, n_nodes] M = parameter_migration_matrices[e];
                W = M' * W * M;
            }

            real t_start;
            real t_end;
            if (e == 1) {
                t_start = 0.0;
                t_end = cumulative_times[1];
            } else if (e == n_events + 1) {
                t_start = cumulative_times[n_events];
                // T_max is the tail length PAST the root, not an absolute cap.
                // It used to be `t_end = T_max`, which silently produced an
                // INVERTED final epoch (t_start > t_end) whenever the root ran
                // deeper than T_max.  On shallow graphs that never bound; on
                // (GBR,IBS,YRI) the SNP term needs t_root of order 1e3 and the
                // (t,Ne) ridge let the optimiser walk to 1e6, past the wall,
                // where int_p_L integrates backwards -- two topologies died on
                // non-finite gradients and the rest fitted garbage.  Anchoring
                // the tail to the root removes the wall entirely and is what
                // T_max was always meant to express: "integrate the deepest
                // epoch far enough for the survival to decay".  For a shallow
                // graph this changes the final epoch by t_root/T_max, i.e. by
                // 0.2% at t_root = 200 -- numerically a no-op.
                t_end = cumulative_times[n_events] + T_max;
            } else {
                t_start = cumulative_times[e - 1];
                t_end = cumulative_times[e];
            }

            real duration = t_end - t_start;

            for (a in 1:n_nodes) {
                real w_diag = W[a, a];

                if (w_diag > 1e-20) {
                    for (b in 1:n_bins) {
                        ibd_fraction[b][li, lj] += w_diag *
                            int_p_L(Ne[a], t_start, t_end,
                                    bin_length[b, 1], bin_length[b, 2]);
                        ibd_number[b][li, lj] += w_diag *
                            int_p_LL(Ne[a], t_start, t_end,
                                     bin_length[b, 1], bin_length[b, 2]);
                    }

                    W[a, a] = w_diag * exp(-duration / Ne[a]);
                }
            }
        }

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
    admixture_fractions ~ beta(1.0, 1.0);

    // Non-centered hierarchical log-normal for Ne: Ne_raw carries the per-node
    // variation, so there is no funnel between Ne and sigma_log.  Each Ne has
    // median pinned at 15000; the marginal prior is mean ~16200, SD ~7400
    // (a log-normal's mean sits above its median).
    mu_log    ~ normal(log(15000), 0.53);  // log-scale location; median Ne = 15000
    // 0.53 = sqrt(trigamma(4)), the log-scale sd of the gamma(4, 4/15000) prior
    // that Nfixed puts on effective_N -- so Nsmooth and Nfixed now express the
    // SAME prior belief about the LEVEL of Ne and differ only in whether it may
    // vary across the graph.  At the old 0.25 the two were not comparable: +-28%
    // at 1 sigma made Ne = 3,300 a ~6-sigma excursion, so SNP-only-Nsmooth on the
    // AFR-EUR two-leaf test reported the prior (Ne = 10,100, g ratio 0.35) while
    // SNP-only-Nfixed reached the data (Ne = 3,323, g ratio 1.00).
    sigma_log ~ normal(0, 0.3);            // half-normal (sigma_log >= 0); per-node spread
    tau       ~ normal(0, 0.3);            // half-normal; log-Ne RW step scale
    Ne_raw    ~ std_normal();

    for (i in 1:n_leaves) {
        for (j in i:n_leaves) {
            for (b in 1:n_bins) {
                if (ibd_hat[b][i, j] > 0) {
                    target += normal_lpdf(ibd_hat[b][i,j] | ibd_fraction[b][i, j], ibd_se[b][i, j]);
                } else {
                    if (i == j) {
                        target += -ibd_number[b][i,j] * cm * n_samples[i] * (n_samples[i] - 1) / 2.0;
                    } else {
                        target += -ibd_number[b][i,j] * cm * n_samples[i] * n_samples[j];
                    }
                }
            }
        }
    }
}

generated quantities {
    real lp_ibd;
    real chi2_ibd;
    int n_ibd_obs;

    // IBD log-likelihood, split out so the composite can be decomposed.  Uses the
    // full normal_lpdf (the model block's `~` drops the -0.5*log(2*pi) constants,
    // which cancel within a component but NOT between two components of different
    // size -- and comparing the two components is the whole point here).
    lp_ibd = 0;
    chi2_ibd = 0;
    n_ibd_obs = 0;
    for (i in 1:n_leaves) {
        for (j in i:n_leaves) {
            for (b in 1:n_bins) {
                if (ibd_hat[b][i, j] > 0) {
                    lp_ibd += normal_lpdf(ibd_hat[b][i,j] | ibd_fraction[b][i,j], ibd_se[b][i,j]);
                    chi2_ibd += square((ibd_hat[b][i,j] - ibd_fraction[b][i,j]) / ibd_se[b][i,j]);
                    n_ibd_obs += 1;
                } else if (i == j) {
                    lp_ibd += -ibd_number[b][i,j] * cm * n_samples[i] * (n_samples[i] - 1) / 2.0;
                } else {
                    lp_ibd += -ibd_number[b][i,j] * cm * n_samples[i] * n_samples[j];
                }
            }
        }
    }
}
