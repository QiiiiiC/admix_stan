data {
  int<lower=0> N_obs;      // number of observations
  real<lower=0> N_pop;       // number of populations
  array[N_obs] real d_mean;
  array[N_obs] real<lower = 0> d_var;
  real<lower = 0> ancestral_T;      // the true ancestral time
  array[N_obs,2] int group;
  array[N_obs] real<lower = 0> adjusted_factor;
}

parameters {
    // complete 7 weights
    real<lower = 0> c1;
    real<lower = 0> c2;
    real<lower = 0> c3;
    real<lower = 0> c4;
    real<lower = 0> c5;
    real<lower = 0> c6;
    real<lower = 0> c7;
    real<lower = 0, upper = 1> f;
}

transformed parameters {
    // real covariance matrix

   real a_11 = c1 + c2 + c5;
   real a_12 = f * c2 + c5;
   real a_13 = c5;
   real a_14 = 0;
   real a_22 = c7 + pow(f,2) * c2 + pow(1-f,2) * c4 + c5;
   real a_23 = (1-f) * c4 + c5;
   real a_24 = 0;
   real a_33 = c3 + c4 + c5;
   real a_34 = 0;
   real a_44 = c6;


   real W_11 = a_11 - 2/N_pop * (a_11 + a_12 + a_13 + a_14) + 1/pow(N_pop,2) * (a_11 + a_22 + a_33 + a_44 + 2 * (a_12 + a_13 + a_14 + a_23 + a_34 + a_24));
   real W_12 = a_12 - 1/N_pop * (a_11 + a_12 + a_13 + a_14) - 1/N_pop * (a_12 + a_22 + a_23 + a_24) + 1/pow(N_pop,2) * (a_11 + a_22 + a_33 + a_44 + 2 * (a_12 + a_13 + a_14 + a_23 + a_34 + a_24));
   real W_13 = a_13 - 1/N_pop * (a_11 + a_12 + a_13 + a_14) - 1/N_pop * (a_13 + a_23 + a_33 + a_34) + 1/pow(N_pop,2) * (a_11 + a_22 + a_33 + a_44 + 2 * (a_12 + a_13 + a_14 + a_23 + a_34 + a_24));
   real W_14 = a_14 - 1/N_pop * (a_11 + a_12 + a_13 + a_14) - 1/N_pop * (a_14 + a_24 + a_34 + a_44) + 1/pow(N_pop,2) * (a_11 + a_22 + a_33 + a_44 + 2 * (a_12 + a_13 + a_14 + a_23 + a_34 + a_24));
   real W_22 = a_22 - 2/N_pop * (a_12 + a_22 + a_23 + a_24) + 1/pow(N_pop,2) * (a_11 + a_22 + a_33 + a_44 + 2 * (a_12 + a_13 + a_14 + a_23 + a_34 + a_24));
   real W_23 = a_23 - 1/N_pop * (a_12 + a_22 + a_23 + a_24) - 1/N_pop * (a_13 + a_23 + a_33 + a_34) + 1/pow(N_pop,2) * (a_11 + a_22 + a_33 + a_44 + 2 * (a_12 + a_13 + a_14 + a_23 + a_34 + a_24));
   real W_24 = a_24 - 1/N_pop * (a_12 + a_22 + a_23 + a_24) - 1/N_pop * (a_14 + a_24 + a_34 + a_44) + 1/pow(N_pop,2) * (a_11 + a_22 + a_33 + a_44 + 2 * (a_12 + a_13 + a_14 + a_23 + a_34 + a_24));
   real W_33 = a_33 - 2/N_pop * (a_13 + a_23 + a_33 + a_34) + 1/pow(N_pop,2) * (a_11 + a_22 + a_33 + a_44 + 2 * (a_12 + a_13 + a_14 + a_23 + a_34 + a_24));
   real W_34 = a_34 - 1/N_pop * (a_13 + a_23 + a_33 + a_34) - 1/N_pop * (a_14 + a_24 + a_34 + a_44)+ 1/pow(N_pop,2) * (a_11 + a_22 + a_33 + a_44 + 2 * (a_12 + a_13 + a_14 + a_23 + a_34 + a_24));
   real W_44 = a_44 - 2/N_pop * (a_14 + a_24 + a_34 + a_44) + 1/pow(N_pop,2) * (a_11 + a_22 + a_33 + a_44 + 2 * (a_12 + a_13 + a_14 + a_23 + a_34 + a_24));

   array[10] real W_theoretical;
    W_theoretical[1] = W_11;  // (0,0)
    W_theoretical[2] = W_12;  // (0,1)
    W_theoretical[3] = W_13;  // (0,2)
    W_theoretical[4] = W_14;  // (0,3)
    W_theoretical[5] = W_22;  // (1,1)
    W_theoretical[6] = W_23;  // (1,2)
    W_theoretical[7] = W_24;  // (1,3)
    W_theoretical[8] = W_33;  // (2,2)
    W_theoretical[9] = W_34;  // (2,3)
    W_theoretical[10] = W_44; // (3,3)

}


model {
    c1 ~ exponential(100);  
    c2 ~ exponential(100);
    c3 ~ exponential(100);
    c4 ~ exponential(100);
    c5 ~ exponential(100);
    c6 ~ exponential(100);
    c7 ~ exponential(100);

    f ~ beta(2,2);
        for (i in 1:N_obs) {
        target += normal_lpdf(d_mean[i] | W_theoretical[i]*100, d_var[i]);
    }

}