functions {
    // Returns the analytical integral of p(l in [u,v] | g) * p(g | N, t1, t2)
    // over g in [t1, t2], where coalescence rate is 1/N (haploid).
    // p(g | same pop) = (1/N) * exp(-(g - t1)/N).
    // length measured in centiMorgans.
    real int_p_L(real N, real t1, real t2, real u, real v) {
        real cu = -(1.0/N + u/50.0);
        real cv = -(1.0/N + v/50.0);
        real k1 = u/(50.0*N) * (exp(cu*t2)*(cu*t2 - 1.0)/pow(cu, 2)
                               - exp(cu*t1)*(cu*t1 - 1.0)/pow(cu, 2));
        real k2 = v/(50.0*N) * (exp(cv*t2)*(cv*t2 - 1.0)/pow(cv, 2)
                               - exp(cv*t1)*(cv*t1 - 1.0)/pow(cv, 2));
        real k3 = 1.0/N * (exp(cu*t2)/cu - exp(cu*t1)/cu);
        real k4 = 1.0/N * (exp(cv*t2)/cv - exp(cv*t1)/cv);

        return exp(t1/N) * (k1 - k2 + k3 - k4);
    }

    vector vc_int_p_L(vector N_vec, real t1, real t2, real u, real v) {
        int len = rows(N_vec);
        vector[len] result;
        for (i in 1:len) {
            result[i] = int_p_L(N_vec[i], t1, t2, u, v);
        }
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
    // admixture_map[i] = [event_index, src, tgt1, tgt2]

    // FIX (index bounds): constrain index arrays so Stan catches bad values
    // at data-load time rather than giving a cryptic runtime error.
    array[n_admixture] int<lower=1, upper=n_events+1> admixture_indices;
    array[n_events - n_admixture] int<lower=1, upper=n_events> fixed_indices;
    array[n_events - n_admixture] int<lower=1, upper=n_events+1> fixed_indices_shifted;

    // Fixed parameter block
    vector<lower=0>[n_nodes] effective_N;   // haploid effective size
    int<lower=0> n_bins;
    array[n_bins, 2] real<lower=0> bin_length;

    // FIX (T_max): user-supplied upper bound for the final open epoch,
    // replacing the hard-coded 1e6. Set to e.g. max(effective_N) * 100
    // or the oldest plausible TMRCA in your system.
    real<lower=0> T_max;

    // Observation block
    array[n_bins] matrix<lower=0>[n_leaves, n_leaves] ibd_hat;
    array[n_bins] matrix<lower=0>[n_leaves, n_leaves] ibd_se;
}

parameters {
    // Incremental epoch durations — must be positive.
    vector<lower=0>[n_events] times;
    array[n_admixture] real<lower=0, upper=1> admixture_fractions;
}

transformed parameters {
    vector[n_events] cumulative_times = cumulative_sum(times);
    matrix[n_nodes, n_nodes] I = diag_matrix(rep_vector(1.0, n_nodes));

    // 1. Build parameter_migration_matrices (identity + fixed + admixture)
    array[n_events + 1] matrix[n_nodes, n_nodes] parameter_migration_matrices;
    parameter_migration_matrices[1] = I;

    // FIX (index bounds): indices are now declared with bounds above,
    // so any out-of-range value will be caught at data loading.
    for (k in 1:(n_events - n_admixture)) {
        parameter_migration_matrices[fixed_indices_shifted[k]] =
            migration_matrices[fixed_indices[k]];
    }

    for (i in 1:n_admixture) {
        int src  = admixture_map[i, 2];
        int tgt1 = admixture_map[i, 3];
        int tgt2 = admixture_map[i, 4];

        matrix[n_nodes, n_nodes] admix_mat = I;
        admix_mat[tgt1, src] =  admixture_fractions[i];
        admix_mat[tgt2, src] =  1.0 - admixture_fractions[i];
        admix_mat[src,  src] =  0.0;

        parameter_migration_matrices[admixture_indices[i]] = admix_mat;
    }

    // 2. Initialise IBD accumulator and lineage-tracking matrices
    array[n_bins] matrix[n_nodes, n_nodes] ibd_fraction =
        rep_array(rep_matrix(0.0, n_nodes, n_nodes), n_bins);

    // current_location[i, j]: probability that a lineage drawn from node i
    // ends up in node j during the current epoch.
    matrix[n_nodes, n_nodes] current_location = I;

    // past_no_coal_prob[i, j]: probability that lineages from nodes i and j
    // have NOT coalesced in any epoch before the current one.
    // At time 0 no pair has had any chance to coalesce → all ones.
    matrix[n_nodes, n_nodes] past_no_coal_prob = rep_matrix(1.0, n_nodes, n_nodes);

    // 3. Loop over epochs
    for (e in 1:(n_events + 1)) {
        if (e > 1) {
            // Propagate lineage locations through the migration/admixture event
            // that separates epoch e-1 from epoch e.
            current_location = current_location * parameter_migration_matrices[e];

            // FIX (factor-of-2): survival probability over epoch e-1 uses
            // haploid coalescence rate 1/N, consistent with int_p_L.
            // Previously this was 1/(2N) (diploid), causing a 2x mismatch.
            vector[n_nodes] lambda = exp(-times[e-1] ./ effective_N);
            matrix[n_nodes, n_nodes] lambda_D = diag_matrix(lambda);

            // Update no-coalescence probability: multiply in the survival
            // probability for epoch e-1, weighted by current lineage locations.
            // quad_form(A, B) = B' * A * B, so with B = current_location':
            //   result[i,j] = sum_k loc[i,k] * lambda[k] * loc[j,k]  ✓
            past_no_coal_prob = past_no_coal_prob .*
                quad_form(lambda_D, current_location');
        }

        // Accumulate IBD contribution from epoch e
        for (b in 1:n_bins) {
            vector[n_nodes] ibd;

            if (e == 1) {
                // First epoch: starts at generation 0
                ibd = vc_int_p_L(effective_N, 0.0, cumulative_times[1],
                                 bin_length[b, 1], bin_length[b, 2]);
            } else if (e == (n_events + 1)) {
                // FIX (T_max): final open epoch now uses data-supplied T_max
                // instead of the hard-coded 1e6, which was fragile for large N.
                ibd = vc_int_p_L(effective_N, cumulative_times[n_events], T_max,
                                 bin_length[b, 1], bin_length[b, 2]);
            } else {
                ibd = vc_int_p_L(effective_N, cumulative_times[e-1],
                                 cumulative_times[e],
                                 bin_length[b, 1], bin_length[b, 2]);
            }

            // Weight epoch IBD by the probability of not having coalesced
            // before this epoch (past_no_coal_prob already updated above).
            ibd_fraction[b] += quad_form(diag_matrix(ibd), current_location') .*
                               past_no_coal_prob;
        }
    }
}

model {
    // FIX (times prior): exponential(0.01) gives mean 100 generations per
    // increment, which compounds badly over many events. A half-normal or
    // exponential with a tighter rate is more appropriate for demographic
    // inference where epoch durations are rarely >1000 generations each.
    // Adjust the scale to match your expected tree height.
    times ~ exponential(0.001);   // mean 1000 gen per increment; tune as needed

    admixture_fractions ~ uniform(0, 1);

    // Likelihood: upper triangle (including diagonal) only, since IBD
    // sharing is symmetric and each pair should be counted once.
    for (i in 1:n_leaves) {
        for (j in i:n_leaves) {
            for (b in 1:n_bins) {
                ibd_hat[b][i, j] ~ normal(ibd_fraction[b][i, j], ibd_se[b][i, j]);
            }
        }
    }
}
