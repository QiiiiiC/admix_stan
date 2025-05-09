data {
  int<lower=0> N_obs;      // number of observations
  int<lower=0> N_pop;       // number of populations
  vector<lower=0>[N_obs] u;     // starting interval
  vector<lower=0>[N_obs] v;     // ending interval
  array[N_obs] int mean_fraction;     // observed mean fraction
  array[N_obs] int mean_number;    // observed mean number of sharing
  array[N_obs,2] int group;     //number of different groups
  real length;      //length of the genome considered measured in cM
  array[N_obs] real<lower = 0> boot_var;     //bootstrap variance measured from data
  real<lower = 0> outgroup_N;       // the N_e of the outgroup
  array[N_obs] real<lower=0, upper = 1> p_0;     //the empirical number of zeros
  real<lower = 0> ancestral_N;      // the N_e of the ancestral population
  real<lower = 0> ancestral_T;      // the true ancestral time
  real<lower = 1> T2;       // the merge time
}

parameters {
    vector<lower=1, upper=100000>[4] N;     // effective population sizes
    real<lower=1, upper = T2> T1;     //admixture times 
    real<lower=0, upper = 1> fraction;      //admixture fraction
}

transformed parameters {
   real c1 = T1/N[1];
   real c2 = (T2-T1)/N[1];
   real c3 = T1/N[3];
   real c4 = (T2-T1)/N[3];
   real c5 = (ancestral_T-T2)/N[4];
}


model {
    for (i in 1:4) {
      N[i] ~ uniform(1,100000);     //prior such that mean is 5000 and variance is 2000
    }
    T1 ~ uniform(1,T2);
    fraction ~ uniform(0,1);
    real c6 = ancestral_T/outgroup_N;

    for (i in 1:N_obs) {
        // The sharing between any population and the outergroup is ignored as there's no parameter involved.
      // This is the sharing between population A and A
      real lambda;
      if (group[i][1]==0 && group[i][2]==0){
        target += normal_lpdf(c1+c2+c5 | lambda, boot_var[i]);
      }
      // This is the sharing between population A and Admix
      else if (group[i][1]==0 && group[i][2]==1){
        lambda = (fraction*int_p_L_divide_l(N[1],T1,ancestral_T,u[i],v[i]) + int_p_L(N[6],ancestral_T,10000000,u[i], v[i]))*100;
        target += normal_lpdf(mean_fraction[i] | lambda, boot_var[i]);
      }
      // This is the sharing between population A and B
      else if (group[i][1]==0 && group[i][2]==2){
        lambda = (int_p_L(N[4],T2, ancestral_T,u[i], v[i]) + int_p_L(ancestral_N, ancestral_T, 10000000, u[i], v[i]))*100;
        target += normal_lpdf(mean_fraction[i] | lambda, boot_var[i]);
      }

      // This is the sharing between population Admix and Admix
      else if (group[i][1]==1 && group[i][2]==1){
        lambda = (int_p_L(N[2],0,T1,u[i],v[i]) + pow(fraction,2) * int_p_L(N[1],T1,ancestral_T,u[i],v[i]) + pow(1-fraction,2)*(int_p_L(N[3],T1,T2,u[i],v[i])+int_p_L(N[4],T2,ancestral_T,u[i],v[i])) + int_p_L(ancestral_N,ancestral_T,100000000,u[i], v[i]))*100;
        target += normal_lpdf(mean_fraction[i] | lambda, boot_var[i]);
      }
      // This is the sharing number between population Admix and B
      else if (group[i][1]==1 && group[i][2]==2){
        lambda = ((1-fraction) * (int_p_L(N[3],T1,T2,u[i],v[i]) + int_p_L(N[4],T2,ancestral_T,u[i],v[i])) + int_p_L(ancestral_N,ancestral_T,100000000,u[i], v[i]))*100;
        target += normal_lpdf(mean_fraction[i] | lambda, boot_var[i]);
      }

      // This is the sharing number between population B and B
      else if (group[i][1]==2 && group[i][2]==2){
        lambda = (int_p_L_divide_l(N[3],0,T2,u[i],v[i]) + int_p_L_divide_l(N[4],T2,ancestral_T,u[i], v[i])+int_p_L(ancestral_N,ancestral_T,10000000,u[i], v[i]))*100;
        target += normal_lpdf(mean_fraction[i] | lambda, boot_var[i]);
      }
    }

}