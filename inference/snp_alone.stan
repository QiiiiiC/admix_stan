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
    vector<lower=1, upper=100000>[5] N;     // effective population sizes
    real<lower=0, upper = 1> fraction;      //admixture fraction
    real <lower = 1, upper = ancestral_T> T2;
    real<lower=1, upper = T2> T1;     //admixture times
}

transformed parameters {
   real c1 = T1/N[1];
   real c2 = (T2-T1)/N[1];
   real c3 = T1/N[3];
   real c4 = (T2-T1)/N[3];
   real c5 = (ancestral_T-T2)/N[4];
   real c6 = ancestral_T/N[5];
   real c7 = T1/N[2];
   real a_11 = c1 + c2 + c5;
   real a_12 = fraction * c2 + c5;
   real a_13 = c5;
   real a_14 = 0;
   real a_22 = c7 + pow(fraction,2) * c2 + pow(1-fraction,2) * c4 + c5;
   real a_23 = (1-fraction)*c4 + c5;
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
}


model {
    for (i in 1:4) {
      N[i] ~ uniform(1,100000);     //prior such that mean is 5000 and variance is 2000
    }
    T2 ~ uniform(1,ancestral_T);
    T1 ~ uniform(1,T2);

    fraction ~ uniform(0,1);

    for (i in 1:N_obs) {
        // The sharing between any population and the outergroup is ignored as there's no parameter involved.
      // This is the sharing between population A and A
      if (group[i][1]==0 && group[i][2]==0){
        target += normal_lpdf(d_mean[i]| W_11*100*adjusted_factor[1], d_var[i]);
      }
      // This is the sharing between population A and Admix
      else if (group[i][1]==0 && group[i][2]==1){
        target += normal_lpdf(d_mean[i]| W_12*100*adjusted_factor[1], d_var[i]);
      }
      // This is the sharing between population A and B
      else if (group[i][1]==0 && group[i][2]==2){
        target += normal_lpdf(d_mean[i]| W_13*100*adjusted_factor[1], d_var[i]);
      }

      // This is the sharing between population Admix and Admix
      else if (group[i][1]==0 && group[i][2]==3){
        target += normal_lpdf(d_mean[i]| W_14*100*adjusted_factor[1], d_var[i]);
      }
      // This is the sharing number between population Admix and B
      else if (group[i][1]==1 && group[i][2]==1){
        target += normal_lpdf(d_mean[i]| W_22*100*adjusted_factor[1], d_var[i]);
      }

      // This is the sharing number between population B and B
      else if (group[i][1]==1 && group[i][2]==2){
        target += normal_lpdf(d_mean[i]| W_23*100*adjusted_factor[1], d_var[i]);
      }

      else if (group[i][1]==1 && group[i][2]==3){
        target += normal_lpdf(d_mean[i]| W_24*100*adjusted_factor[1], d_var[i]);
      }
      else if (group[i][1]==2 && group[i][2]==2){
        target += normal_lpdf(d_mean[i]| W_33*100*adjusted_factor[1], d_var[i]);
      }
      else if (group[i][1]==2 && group[i][2]==3){
        target += normal_lpdf(d_mean[i]| W_34*100*adjusted_factor[1], d_var[i]);
      }
      else if (group[i][1]==3 && group[i][2]==3){
        target += normal_lpdf(d_mean[i]| W_44*100*adjusted_factor[1], d_var[i]);
      }
    }

}