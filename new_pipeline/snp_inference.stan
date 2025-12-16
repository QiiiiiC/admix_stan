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
  vector<lower = 0>[n_nodes] effective_N;

    // Observation block
  matrix[n_leaves, n_leaves] w_hat;
  matrix<lower=0>[n_leaves, n_leaves] w_se;
}

parameters{
    // vector<lower = 0>[n_nodes] effective_N;   // effective population sizes, assume constant for each nodes
    vector<lower = 0>[n_events] times;          // times between each events
    array[n_admixture] real<lower = 0, upper = 1> admixture_fractions; // admixture fractions for each admixture event(default to the first population)
}

transformed parameters{
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

    matrix[n_nodes, n_nodes] V = rep_matrix(0, n_nodes, n_nodes);
    matrix[n_nodes, n_nodes] current_location = I; 
    

    for (e in 1:n_events) {
        current_location = current_location * parameter_migration_matrices[e];
        vector[n_nodes] drift_factors = times[e] ./ (2.0 * effective_N);
        matrix[n_nodes, n_nodes] D = diag_matrix(drift_factors);
        V += quad_form(D, current_location');
    }   // Compute the theoretical mean for covariance matrix


    matrix[n_nodes, n_nodes] W;
    {
        vector[n_nodes] row_means;
        real grand_mean;

        // Calculate Row Means
        for (i in 1:n_nodes) {
            row_means[i] = mean(V[i]);
        }
        
        // Calculate Grand Mean
        grand_mean = mean(row_means);

        W = V - rep_matrix(row_means, n_nodes) 
              - rep_matrix(row_means', n_nodes) 
              + grand_mean;
    }
}

model {
    times ~ exponential(1.0/100);  
    admixture_fractions ~ uniform(0,1);   

    for (i in 1:n_leaves) {
        for (j in i:n_leaves) {
            w_hat[i, j] ~ normal(W[i, j], w_se[i, j]);
        }
    }
}