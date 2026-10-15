// Mixed (IBD + SNP) model, smoothed Ne, Poisson IBD counts -- per-branch formulation.
//
// Same likelihood and the same priors as methods/stan/mixed_model_Nsmooth_poisson.stan, but
// the graph enters as a TOPOLOGY (a partial order of events), not as an event order:
//
//   * `times[e]` is the wait of event e after the latest of the events it must come
//     after; event heights follow from the topology.  Unrelated events are not ordered.
//   * every branch contributes once, over its own lifespan [start, end): its survival
//     exp(-(end - start) / Ne) and its segment density per bin.  No epochs, no sorting.
//   * one sweep over the events (in any fixed listing that respects the partial order)
//     tracks T[x, y][i, j] = P(lineage of leaf i in branch x, lineage of leaf j in branch y,
//     not coalesced), for all leaf pairs at once.
//   * expected IBD = (segment density by bin and branch) x (mass arriving in each branch, by
//     leaf pair): one matrix product.
//   * SNP drift = sum over branches of p_i(x) p_j(x) length / Ne, with p_i(x) the
//     probability that leaf i's lineage passes through x.
data {
    int<lower=1> n_leaves;
    int<lower=1> n_nodes;                     // branches; leaves first, a branch's parent has a higher index
    int<lower=1> n_events;
    int<lower=0> n_admixture;

    // ---- Events, in any listing that respects "must come after" ----
    array[n_events] int<lower=0, upper=n_admixture> event_admix;   // admixture index; 0 = merge
    array[n_events, 3] int<lower=1, upper=n_nodes> event_nodes;    // merge: child, child, parent
                                                                   // admixture: child, source (f), source (1 - f)
    array[n_events] int<lower=0, upper=2> n_pred;                  // events this one must come after
    array[n_events, 2] int<lower=0, upper=n_events> pred;
    array[n_events] int<lower=0, upper=n_nodes> n_open;            // other branches open when it is processed
    array[n_events, n_nodes] int<lower=0, upper=n_nodes> open;
    array[n_admixture] int<lower=1, upper=n_events> admix_event;

    // ---- Branches ----
    array[n_nodes] int<lower=0, upper=n_events> start_event;       // creating event; 0 = leaf (starts at 0)
    array[n_nodes] int<lower=0, upper=n_events> end_event;         // consuming event; 0 = root
    int<lower=1, upper=n_nodes> root;
    array[n_nodes] int<lower=0> ne_parent;                         // branch it flows into; 0 = root / admixed
    array[n_nodes] int<lower=0> ne_admix_idx;                      // admixture index if the branch is admixed

    array[n_leaves] int<lower=0> n_samples;

    // ---- IBD ----
    int<lower=0> n_bins;
    array[n_bins, 2] real<lower=0> bin_length;
    real<lower=0> T_max;                                           // length of the root branch
    real<lower=0> cm;
    array[n_bins, n_leaves, n_leaves] int<lower=0> ibd_count;

    // ---- SNP ----
    matrix[n_leaves, n_leaves] w_hat;
    matrix<lower=0>[n_leaves, n_leaves] w_se;
}

transformed data {
    // Leaf pairs i <= j, and the distinct bin edges (adjacent bins share one).
    int n_pairs = (n_leaves * (n_leaves + 1)) %/% 2;
    array[n_pairs] int pair_i;
    array[n_pairs] int pair_j;
    {
        int p = 0;
        for (i in 1:n_leaves) for (j in i:n_leaves) { p += 1; pair_i[p] = i; pair_j[p] = j; }
    }
    int n_edges = 0;
    array[n_bins] int lo;                     // edge index of each bin's lower and upper bound
    array[n_bins] int hi;
    vector[2 * n_bins] edge_all;
    for (b in 1:n_bins) {
        for (s in 1:2) {
            int found = 0;
            for (k in 1:n_edges) if (abs(edge_all[k] - bin_length[b, s]) < 1e-12) found = k;
            if (found == 0) { n_edges += 1; edge_all[n_edges] = bin_length[b, s]; found = n_edges; }
            if (s == 1) lo[b] = found; else hi[b] = found;
        }
    }
    vector[n_edges] u50 = edge_all[1:n_edges] / 50.0;          // recombination rate per generation
    vector[n_edges] root_tail = exp(-u50 * T_max);
    // IBD counts and exposure as one vector over (bin, pair), bins fastest.
    array[n_bins * n_pairs] int count_flat;
    vector[n_bins * n_pairs] exposure;
    vector[n_pairs] w_obs;
    vector[n_pairs] w_obs_se;
    for (p in 1:n_pairs) {
        int i = pair_i[p];
        int j = pair_j[p];
        real pairs = (i == j) ? n_samples[i] * (n_samples[i] - 1) / 2.0 : n_samples[i] * 1.0 * n_samples[j];
        for (b in 1:n_bins) {
            count_flat[(p - 1) * n_bins + b] = ibd_count[b, i, j];
            exposure[(p - 1) * n_bins + b] = cm * pairs;
        }
        w_obs[p] = w_hat[i, j];
        w_obs_se[p] = w_se[i, j];
    }
}

parameters {
    vector<lower=1>[n_events] times;          // wait after the latest predecessor event
    array[n_admixture] real<lower=0, upper=1> admixture_fractions;
    real mu_log;
    real<lower=0> sigma_log;
    real<lower=0> tau;
    vector[n_nodes] Ne_raw;
}

transformed parameters {
    // ---- Event heights from the partial order (named as in the event-order models) ----
    vector[n_events] cumulative_times;
    for (e in 1:n_events) {
        real base = 0;
        for (k in 1:n_pred[e]) base = fmax(base, cumulative_times[pred[e, k]]);
        cumulative_times[e] = base + times[e];
    }
    vector[n_nodes] node_t;                   // branch start
    vector[n_nodes] node_end;                 // branch end
    for (a in 1:n_nodes) {
        node_t[a] = start_event[a] == 0 ? 0.0 : cumulative_times[start_event[a]];
        node_end[a] = end_event[a] == 0 ? node_t[a] + T_max : cumulative_times[end_event[a]];
    }

    // ---- Smoothed Ne: log-Ne random walk along the flow, as in mixed_model_Nsmooth_poisson ----
    // Each branch is centred on the branch it flows into; an admixed branch on the
    // fraction-weighted mean of its two sources.  Step variance = tau^2 * (own length) / t_ref.
    real t_ref = 100.0;
    vector[n_nodes] log_Ne;
    for (ii in 1:n_nodes) {
        int a = n_nodes - ii + 1;
        real dt = node_end[a] - node_t[a];
        if (ne_parent[a] == 0 && ne_admix_idx[a] == 0) {
            log_Ne[a] = mu_log + sigma_log * Ne_raw[a];
        } else if (ne_admix_idx[a] > 0) {
            int e = admix_event[ne_admix_idx[a]];
            real f = admixture_fractions[ne_admix_idx[a]];
            log_Ne[a] = f * log_Ne[event_nodes[e, 2]] + (1.0 - f) * log_Ne[event_nodes[e, 3]]
                        + tau * sqrt(dt / t_ref + 1e-9) * Ne_raw[a];
        } else {
            log_Ne[a] = log_Ne[ne_parent[a]] + tau * sqrt(dt / t_ref + 1e-9) * Ne_raw[a];
        }
    }
    vector<lower=0>[n_nodes] Ne = exp(log_Ne);

    array[n_bins] matrix[n_leaves, n_leaves] ibd_number;
    matrix[n_leaves, n_leaves] W_centered;
    {
        // ---- 1. Per-branch quantities: each uses only its own start, end and Ne ----
        // exp(-u * height / 50) once per event and bin edge; every branch starts and ends at an event.
        array[n_events] vector[n_edges] decay;
        for (e in 1:n_events) decay[e] = exp(-u50 * cumulative_times[e]);
        vector[n_nodes] surv = exp(-(node_end - node_t) ./ Ne);
        matrix[n_bins, n_nodes] dens;         // segments per pair per cM, by bin and branch
        for (a in 1:n_nodes) {
            // Antiderivative of the segment density at each edge (closed form of int_p_LL).
            vector[n_edges] c = u50 + 1.0 / Ne[a];
            vector[n_edges] g;
            if (start_event[a] == 0) g = c * node_t[a] + 1.0;
            else g = decay[start_event[a]] .* (c * node_t[a] + 1.0);
            if (end_event[a] == 0) {
                if (start_event[a] == 0) g -= surv[a] * (root_tail .* (c * node_end[a] + 1.0));
                else g -= surv[a] * (decay[start_event[a]] .* root_tail .* (c * node_end[a] + 1.0));
            } else {
                g -= surv[a] * (decay[end_event[a]] .* (c * node_end[a] + 1.0));
            }
            g = g ./ square(c) / (50.0 * Ne[a]);
            dens[:, a] = g[lo] - g[hi];
        }

        // ---- 2. One sweep: uncoalesced joint mass for all leaf pairs; route probabilities ----
        array[n_nodes, n_nodes] matrix[n_leaves, n_leaves] T;
        matrix[n_nodes, n_pairs] arrive;      // mass entering each branch with both lineages in it
        matrix[n_nodes, n_leaves] route = rep_matrix(0.0, n_nodes, n_leaves);
        for (i in 1:n_leaves) {
            route[i, i] = 1.0;
            for (j in 1:n_leaves) T[i, j] = rep_matrix(0.0, n_leaves, n_leaves);
            for (j in 1:n_leaves) T[i, j][i, j] = 1.0;
        }
        for (e in 1:n_events) {
            int x = event_nodes[e, 1];
            int y = event_nodes[e, 2];
            int z = event_nodes[e, 3];
            if (event_admix[e] == 0) {                       // merge x, y -> z
                T[z, z] = surv[x] * T[x, x] + surv[y] * T[y, y] + T[x, y] + T[y, x];
                for (k in 1:n_open[e]) {
                    int w = open[e, k];
                    T[z, w] = T[x, w] + T[y, w];
                    T[w, z] = T[w, x] + T[w, y];
                }
                route[z] = route[x] + route[y];
            } else {                                         // admixture x -> y (f), z (1 - f)
                real f = admixture_fractions[event_admix[e]];
                matrix[n_leaves, n_leaves] S = surv[x] * T[x, x];
                T[y, y] = square(f) * S;
                T[z, z] = square(1.0 - f) * S;
                T[y, z] = f * (1.0 - f) * S;
                T[z, y] = T[y, z];
                for (k in 1:n_open[e]) {
                    int w = open[e, k];
                    T[y, w] = f * T[x, w];
                    T[z, w] = (1.0 - f) * T[x, w];
                    T[w, y] = f * T[w, x];
                    T[w, z] = (1.0 - f) * T[w, x];
                }
                route[y] = f * route[x];
                route[z] = (1.0 - f) * route[x];
            }
        }
        for (a in 1:n_nodes) for (p in 1:n_pairs) arrive[a, p] = T[a, a][pair_i[p], pair_j[p]];

        // Expected IBD = sum over branches of (segment density) x (mass arriving in the branch).
        matrix[n_bins, n_pairs] number = dens * arrive;
        for (p in 1:n_pairs) {
            for (b in 1:n_bins) {
                ibd_number[b][pair_i[p], pair_j[p]] = number[b, p];
                ibd_number[b][pair_j[p], pair_i[p]] = number[b, p];
            }
        }

        // ---- 3. SNP drift up to the root, then double centring ----
        vector[n_nodes] drift = (node_end - node_t) ./ Ne;
        drift[root] = 0;
        matrix[n_leaves, n_leaves] snp_mean = route' * diag_pre_multiply(drift, route);
        vector[n_leaves] row_means;
        for (i in 1:n_leaves) row_means[i] = mean(snp_mean[i]);
        real grand_mean = mean(row_means);
        W_centered = snp_mean - rep_matrix(row_means, n_leaves) - rep_matrix(row_means', n_leaves) + grand_mean;
    }
}

model {
    times ~ exponential(0.01);
    admixture_fractions ~ beta(1.0, 1.0);
    mu_log    ~ normal(log(15000), 0.53);
    sigma_log ~ normal(0, 0.3);
    tau       ~ normal(0, 0.3);
    Ne_raw    ~ std_normal();

    {
        vector[n_bins * n_pairs] number;      // per (bin, pair), bins fastest
        for (p in 1:n_pairs)
            for (b in 1:n_bins) number[(p - 1) * n_bins + b] = ibd_number[b][pair_i[p], pair_j[p]];
        target += poisson_lpmf(count_flat | number .* exposure);
        vector[n_pairs] w_model;
        for (p in 1:n_pairs) w_model[p] = W_centered[pair_i[p], pair_j[p]];
        w_obs ~ normal(w_model, w_obs_se);
    }
}

generated quantities {
    real lp_ibd = 0;
    real chi2_ibd = 0;
    real lp_snp = 0;
    real chi2_snp = 0;
    for (p in 1:n_pairs) {
        int i = pair_i[p];
        int j = pair_j[p];
        for (b in 1:n_bins) {
            int k = (p - 1) * n_bins + b;
            real lam = fmax(ibd_number[b][i, j] * exposure[k], 1e-12);
            lp_ibd += poisson_lpmf(count_flat[k] | lam);
            chi2_ibd += square(count_flat[k] - lam) / lam;
        }
        lp_snp += normal_lpdf(w_obs[p] | W_centered[i, j], w_obs_se[p]);
        chi2_snp += square((w_obs[p] - W_centered[i, j]) / w_obs_se[p]);
    }
}
