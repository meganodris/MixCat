data {

 int N;                             // N individuals
 vector[N] y;                       // antibody titers
 array[N] int ageG;                 // age group index
 int NAgeG;                         // N age groups
 vector[NAgeG] ageL;                // lower bound of each age group
 vector[NAgeG] ageU;                // upper bound of each age group
 int NGroup;                        // N foi groups (1 = unstratified)
 array[N] int group;                 // foi group index for each individual
 array[NGroup, NAgeG] int n_age;    // individuals per foi group x age group
 int predL;                         // length of titer fit values
 vector[predL] y_fit;               // titer fit values

 vector[4] prior_means;
 vector[4] prior_sds;

}


parameters {

 vector<lower=0>[NGroup] lambda;  // foi per group
 real mu0;                        // mean seroneg
 real <lower=0> mu1;              // mean seropos
 real <lower=0> sd0;              // sd seroneg
 real <lower=0> sd1;              // sd seropos

}


transformed parameters {

  matrix[NAgeG, NGroup] sero;
  array[2] vector[N] pC;
  vector[N] log_lik;


  //--- seroprevalence per age group x foi group ---//
  for(g in 1:NGroup){
    for(a in 1:NAgeG){
      sero[a,g] = 1 - (exp(-lambda[g] * ageL[a]) - exp(-lambda[g] * ageU[a])) /
                      (lambda[g] * (ageU[a] - ageL[a]));
    }
  }

  //--- likelihood calculation ---//
  for(n in 1:N){
    pC[1,n] = log(1 - sero[ageG[n], group[n]]) + normal_lpdf(y[n] | mu0, sd0);
    pC[2,n] = log(sero[ageG[n], group[n]])      + normal_lpdf(y[n] | mu0+mu1, sd1);
    log_lik[n] = log_sum_exp(pC[,n]);
  }
}


model {

 // priors
 lambda ~ exponential(4);
 mu0 ~ normal(prior_means[1], prior_sds[1]);
 mu1 ~ normal(prior_means[2], prior_sds[2]);
 sd0 ~ normal(prior_means[3], prior_sds[3]);
 sd1 ~ normal(prior_means[4], prior_sds[4]);

 // log-likelihood
 target += sum(log_lik);

}


generated quantities {

  array[N] int z;          // latent component: 0 = seroneg, 1 = seropos
  vector[N] y_rep;         // posterior predictive titer
  vector[NGroup] sero_grp; // seroprevalence per foi group (age-weighted)
  real sero_all;           // overall seroprevalence across all groups and age groups
  vector[predL] fitNeg;
  vector[predL] fitPos;
  vector[predL] fitAll;


  for(n in 1:N){
    real prob1 = exp(pC[2,n] - log_lik[n]);

    z[n] = binomial_rng(1, prob1);
    y_rep[n] = z[n] == 0 ? normal_rng(mu0, sd0)
                          : normal_rng(mu0 + mu1, sd1);
  }

  {
    real weighted_sum = 0;
    for(g in 1:NGroup){
      sero_grp[g] = dot_product(to_vector(n_age[g]), col(sero, g)) / sum(n_age[g]);
      weighted_sum += sero_grp[g] * sum(n_age[g]);
    }
    sero_all = weighted_sum / N;
  }

  for(i in 1:predL){
    fitNeg[i] = (1 - sero_all) * exp(normal_lpdf(y_fit[i] | mu0, sd0));
    fitPos[i] = sero_all       * exp(normal_lpdf(y_fit[i] | mu0+mu1, sd1));
    fitAll[i] = fitNeg[i] + fitPos[i];
  }

}
