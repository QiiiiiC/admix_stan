functions {
    real p_L_divide_l(real N, real t1, real t2, real l) {
        real c = -(l / 50 + 1 / N);
        real output = exp(c * t2) * (pow(c * t2 , 2) - 2*c*t2 + 2)/(pow(c,3)) - exp(c * t1) * (pow(c * t1, 2) - 2*c*t1 + 2)/(pow(c,3));
        real r = exp(t1/N)/(2500 * N);
        return r*output;
    }
    real int_p_L_divide_l(real N, real t1, real t2, real u, real v) {
        real c1 = -(v/50 + 1/N);
        real c2 = -(u/50 + 1/N);
        real k1 = (-1/(50*N)) * ((exp(c1*t2)*(c1*t2-1)/pow(c1,2)) - (exp(c1*t1)*(c1*t1-1)/pow(c1,2)));
        real k2 = (1/(50*N)) * ((exp(c2*t2)*(c2*t2-1)/pow(c2,2)) - (exp(c2*t1)*(c2*t1-1)/pow(c2,2)));
        real r = exp(t1/N);
        return r*(k1 + k2);
    }
    real int_p_L_divide_l_limit(real N, real t, real u, real v) {
        real cu = -(1/N + u/50);
        real cv = -(1/N + v/50);
        return -1/(50*N)*exp(t/N) * (exp(cu*t)*(cu*t-1)/pow(cu,2) - exp(cv*t)*(cv*t-1)/pow(cv,2));
    }
    real int_p_L(real N, real t1, real t2, real u, real v) {
        real cu = -(1/N + u/50);
        real cv = -(1/N + v/50);
        real k1 = u/(50*N) * (exp(cu*t2)*(cu*t2-1)/(pow(cu,2)) - exp(cu*t1)*(cu*t1-1)/(pow(cu,2)));
        real k2 = v/(50*N) * (exp(cv*t2)*(cv*t2-1)/(pow(cv,2)) - exp(cv*t1)*(cv*t1-1)/(pow(cv,2)));
        real k3 = 1/N * (exp(cu*t2)/cu - exp(cu*t1)/cu);
        real k4 = 1/N * (exp(cv*t2)/cv - exp(cv*t1)/cv);
        real r = exp(t1/N);
        return r * (k1-k2+k3-k4);
    }
}

data {
  int<lower=0> N_obs;      // number of observations
  vector<lower=0>[N_obs] u;     // starting interval
  vector<lower=0>[N_obs] v;     // ending interval
  array[N_obs] real mean_fraction;     // observed mean fraction
  array[N_obs] real mean_number;    // observed mean number of sharing
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


model {
    for (i in 1:4) {
      N[i] ~ uniform(1,100000);     //prior such that mean is 5000 and variance is 2000
    }
    T1 ~ uniform(1,T2);
    fraction ~ uniform(0,1);

    for (i in 1:N_obs) {
        // The sharing between any population and the outergroup is ignored as there's no parameter involved.
      // This is the sharing between population A and A
      real lambda;
      if (group[i][1]==0 && group[i][2]==0){
        lambda = (int_p_L(N[1],0,T2,u[i],v[i]) + int_p_L(N[4],T2,ancestral_T,u[i],v[i]) + int_p_L(ancestral_N,ancestral_T,10000000,u[i], v[i])) * 1000;
        target += normal_lpdf(mean_fraction[i] | lambda, boot_var[i]);
      }
      // This is the sharing between population A and Admix
      else if (group[i][1]==0 && group[i][2]==1){
        lambda = (fraction*int_p_L_divide_l(N[1],T1,ancestral_T,u[i],v[i]) + int_p_L(N[6],ancestral_T,10000000,u[i], v[i]))*1000;
        target += normal_lpdf(mean_fraction[i] | lambda, boot_var[i]);
      }
      // This is the sharing between population A and B
      else if (group[i][1]==0 && group[i][2]==2){
        lambda = (int_p_L(N[4],T2, ancestral_T,u[i], v[i]) + int_p_L(ancestral_N, ancestral_T, 10000000, u[i], v[i]))*1000;
        target += normal_lpdf(mean_fraction[i] | lambda, boot_var[i]);
      }

      // This is the sharing between population Admix and Admix
      else if (group[i][1]==1 && group[i][2]==1){
        lambda = (int_p_L(N[2],0,T1,u[i],v[i]) + pow(fraction,2) * int_p_L(N[1],T1,ancestral_T,u[i],v[i]) + pow(1-fraction,2)*(int_p_L(N[3],T1,T2,u[i],v[i])+int_p_L(N[4],T2,ancestral_T,u[i],v[i])) + int_p_L(ancestral_N,ancestral_T,100000000,u[i], v[i]))*1000;
        target += normal_lpdf(mean_fraction[i] | lambda, boot_var[i]);
      }
      // This is the sharing number between population Admix and B
      else if (group[i][1]==1 && group[i][2]==2){
        lambda = ((1-fraction) * (int_p_L(N[3],T1,T2,u[i],v[i]) + int_p_L(N[4],T2,ancestral_T,u[i],v[i])) + int_p_L(ancestral_N,ancestral_T,100000000,u[i], v[i]))*1000;
        target += normal_lpdf(mean_fraction[i] | lambda, boot_var[i]);
      }

      // This is the sharing number between population B and B
      else if (group[i][1]==2 && group[i][2]==2){
        lambda = (int_p_L_divide_l(N[3],0,T2,u[i],v[i]) + int_p_L_divide_l(N[4],T2,ancestral_T,u[i], v[i])+int_p_L(ancestral_N,ancestral_T,10000000,u[i], v[i]))*1000;
        target += normal_lpdf(mean_fraction[i] | lambda, boot_var[i]);
      }
    }

}