data {
    // Topology block
  int<lower = 0> n_leaf;       // number of leaf populations
  int<lower = 0> n_event;      // number of total events
  int<lower = 0> n_nodes;     // number of total nodes
  int<lower = 0> n_admixture; // number of admixture events
  array[n_event] matrix[n_nodes, n_nodes] migration_matrices; // migration matrices for each event
  array[n_admixture, 4] int<lower = 0> admixture_map; // for each admixture event, the first item is the index of the event, 
                                                        //the second is the index of the source population, 
                                                        //the third and the fourth indices are the two target populations.

    // Observation block
  int<lower = 0> n_snp_obs;      // number of observations
  array[n_snp_obs] real data_snp_mean;
  array[n_snp_obs] real<lower = 0> data_snp_var;
  real<lower = 0> ancestral_T;      // the true ancestral time
  array[n_snp_obs,2] int group;
  array[n_snp_obs] real<lower = 0> adjusted_factor; // the average of the observed snp frequencies
}

parameters{
    vector<lower = 0>[n_nodes] effective_N;   // effective population sizes, assume constant for each nodes
    vector<lower = 0>[n_event] times;          // times of each event
    array[n_admixture] real<lower = 0, upper = 1> admixture_fractions; // admixture fractions for each admixture event(default to the first population)
}

transformed parameters{
    array[n_admixture] matrix[n_nodes, n_nodes] admixture_matrices; // admixture matrices for each admixture event
    matrix[n_nodes, n_nodes] I = diag_matrix(rep_vector(1.0, n_nodes));
    for (i in 1:n_admixture){
        int event_index = admixture_map[i,1];
        int source_index = admixture_map[i,2];
        int target1_index = admixture_map[i,3];
        int target2_index = admixture_map[i,4];
        admixture_matrices[i] = I;
        admixture_matrices[i][target1_index, source_index] = admixture_fractions[i];
        admixture_matrices[i][target2_index, source_index] = 1 - admixture_fractions[i];
        admixture_matrices[i][source_index, source_index] = 0.0;
    }
    matrix[n_nodes, n_nodes] effective_N_matrix = diag_matrix(effective_N);
}