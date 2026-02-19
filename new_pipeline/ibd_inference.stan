functions {
    real int_p_L(real N, real t1, real t2, real u, real v) {
        real cu = -(1/N + u/50);
        real cv = -(1/N + v/50);
        real k1 = u/(50*N) * (exp(cu*t2)*(cu*t2-1)/(pow(cu,2)) - exp(cu*t1)*(cu*t1-1)/(pow(cu,2)));
        real k2 = v/(50*N) * (exp(cv*t2)*(cv*t2-1)/(pow(cv,2)) - exp(cv*t1)*(cv*t1-1)/(pow(cv,2)));
        real k3 = 1/N * (exp(cu*t2)/cu - exp(cu*t1)/cu);
        real k4 = 1/N * (exp(cv*t2)/cv - exp(cv*t1)/cv);
        return (k1-k2+k3-k4);
    }
    vector vc_int_p_L(vector N_vec, real t1, real t2, real u, real v) {
        int len = rows(N_vec);
        vector[len] result;

        for (i in 1:len) {
            // Calls the scalar function defined above for each element
            result[i] = int_p_L(N_vec[i], t1, t2, u, v);
        }
        return result;
    }
}

data {
    // Topology block
  int<lower = 0> n_leaves;       // number of leaf populations
  int<lower = 0> n_events;      // number of total events, including the default state as the Identity matrix
  int<lower = 0> n_nodes;     // number of total populations
  int<lower = 0> n_admixture; // number of admixture events
  array[n_events] matrix[n_nodes, n_nodes] migration_matrices; // migration matrices for each event
  array[n_admixture, 4] int<lower = 0> admixture_map; // for each admixture event, the first item is the index of the event, 
                                                        //the second is the index of the source population, 
                                                        //the third and the fourth indices are the two target populations.
  array[n_admixture] int<lower = 0> admixture_indices;
  array[n_events - n_admixture] int <lower = 0> fixed_indices;
  array[n_events - n_admixture] int<lower = 0> fixed_indices_shifted;

    // Fixed parameter block
  vector<lower = 0>[n_nodes] effective_N;
  int<lower = 0> n_bins;
  array[n_bins, 2] real<lower = 0> bin_length;

    // Observation block
  array[n_bins] matrix<lower = 0>[n_leaves, n_leaves] ibd_hat;
  array[n_bins] matrix<lower=0>[n_leaves, n_leaves] ibd_se;
}

parameters{
    //vector<lower = 0>[n_nodes] effective_N;   // effective population sizes, assume constant for each nodes
    vector<lower = 0>[n_events] times;          // times between each epoch, the first time is between 0 and the first event.
    array[n_admixture] real<lower = 0, upper = 1> admixture_fractions; // admixture fractions for each admixture event(default to the first population)
}

transformed parameters{
    vector[n_events] cumulative_times = cumulative_sum(times);


    array[n_admixture] matrix[n_nodes, n_nodes] admixture_matrices;
    matrix[n_nodes, n_nodes] I = diag_matrix(rep_vector(1.0, n_nodes));

    for (i in 1:n_admixture){
        int source_index = admixture_map[i,2];
        int target1_index = admixture_map[i,3];
        int target2_index = admixture_map[i,4];
        admixture_matrices[i] = I;
        admixture_matrices[i][target1_index, source_index] = admixture_fractions[i];
        admixture_matrices[i][target2_index, source_index] = 1 - admixture_fractions[i];
        admixture_matrices[i][source_index, source_index] = 0.0;
    }

    array[n_events + 1] matrix[n_nodes, n_nodes] parameter_migration_matrices;
    parameter_migration_matrices[1] = I;
    parameter_migration_matrices[fixed_indices_shifted] = migration_matrices[fixed_indices];
    parameter_migration_matrices[admixture_indices] = admixture_matrices;   // Update the migration matrices with admixture fractions

    array[n_bins] matrix[n_nodes, n_nodes] ibd_fraction = rep_matrix(0, n_nodes, n_nodes);
    matrix[n_nodes, n_nodes] current_location = I; 
    matrix[n_nodes, n_nodes] past_coal_probability = diag_matrix(rep_vector(1.0, n_nodes));
    

    for (e in 1:(n_events+1)) {
        if (e > 1){
            current_location = current_location * parameter_migration_matrices[e];// the i,j entry is the probability of individual drawn from node i ends up in node j between event e and e+1
            vector[n_nodes] lambda = exp(-times[e-1]./ (2*effective_N)) * exp(-cumulative_times[e-1]./ (2* effective_N));
            matrix[n_nodes, n_nodes] lambda_D = diag_matrix(lambda);
            past_no_coal_probability += quad_form(lambda_D, current_location');
        }

        for (b in 1:n_bins){
            if (e == 1){
                vector[n_nodes] ibd = vc_int_p_L(effective_N, 0, cumulative_times[e], bin_length[b,1], bin_length[b,2]);
            }
            else if(e == (n_events + 1)){
                vector[n_nodes] ibd = vc_int_p_L(effective_N, cumulative_times[e-1], 100000000000, bin_length[b,1], bin_length[b,2]);
            }
            else{
                vector[n_nodes] ibd = vc_int_p_L(effective_N, cumulative_times[e-1], cumulative_times[e], bin_length[b,1], bin_length[b,2]);
            }
            ibd_fraction[b] += quad_form(diag_matrix(ibd), current_location') .* past_coal_probability;
        }
        
    }   


    
}

model {
    times ~ exponential(1.0/100);  
    admixture_fractions ~ uniform(0,1);   

    for (i in 1:n_leaves) {
        for (j in i:n_leaves) {
            for (b in 1:n_bins){
                ibd_hat[i, j] ~ normal(ibd_fraction[b][i, j], ibd_se[b][i, j]);
            }
        }
    }
}