data {
  // number of observations
  int<lower=0> N;
  // number of columns in 1st design matrix, including intercept
  int<lower=1> P_1;
  // design matrix for non-terminal model
  matrix[N, P_1] X_1;
  // number of columns in 2nd design matrix, including intercept
  int<lower=1> P_2;
  // design matrix for terminal model w/o non-terminal history
  matrix[N, P_2] X_2;
  // number of columns in 3rd design matrix, including intercept
  int<lower=1> P_3;
  // design matrix for terminal model with non-terminal history
  matrix[N, P_3] X_3;
  // observed non-terminal time
  vector<lower=0>[N] Yr;
  // indicator of event observation for non-terminal event
  int<lower=0,upper=1> dYr[N];
  // observed terminal time
  vector<lower=0>[N] Yt;
  // indicator of event observation for terminal event
  int<lower=0,upper=1> dYt[N];
}

transformed data {
  // Duration of gap between non-terminal and terminal events
  vector[N] YtYrdiff;
  YtYrdiff = Yt - Yr;
}

parameters {
  // vectors of regression parameters
  vector[P_1] beta1;
  vector[P_2] beta2;
  vector[P_3] beta3;

  // shape parameters (the one in exponent of time)
  // alpha > 1 -> hazard increases over time, more clumping
  real<lower=0> alpha1;
  real<lower=0> alpha2;
  real<lower=0> alpha3;
  
  // scale parameters
  // bigger sigma -> slower event occurrence (double check this)
  //real<lower=0> sigma1;
  //real<lower=0> sigma2;
  //real<lower=0> sigma3;
}

model {
  // linear predictors
  vector[N] lp1;
  vector[N] lp2;
  vector[N] lp3;
  lp1 = X_1 * beta1;
  lp2 = X_2 * beta2;
  lp3 = X_3 * beta3;
  
  // no priors -> use Stan defaults
  
  // likelihood
  for (n in 1:N){
    if (dYr[n] == 0 && dYt[n] == 0) {
      
      // type 1: observe neither event
      target += weibull_lccdf(Yr[n] | alpha1, exp(-(lp1[n])/alpha1)) + 
                weibull_lccdf(Yt[n] | alpha2, exp(-(lp2[n])/alpha2));
                
    } else if (dYr[n] == 1 && dYt[n] == 0) {
      
      // type 2: observe non-terminal but terminal censored
      target += weibull_lpdf(Yr[n] | alpha1, exp(-(lp1[n])/alpha1)) +
                weibull_lccdf(Yr[n] | alpha2, exp(-(lp2[n])/alpha2)) + 
                weibull_lccdf(YtYrdiff[n] | alpha3, exp(-(lp3[n])/alpha3));
      
    } else if (dYr[n] == 0 && dYt[n] == 1) {
      
      // type 3: observed terminal with no prior non-terminal
      target += weibull_lccdf(Yr[n] | alpha1, exp(-(lp1[n])/alpha1)) +
                weibull_lpdf(Yt[n] | alpha2, exp(-(lp2[n])/alpha2));
      
    } else if (dYr[n] == 1 && dYt[n] == 1) {
      
      // type 4: both non-terminal and terminal observed
      target += weibull_lpdf(Yr[n] | alpha1, exp(-(lp1[n])/alpha1)) +
                weibull_lccdf(Yr[n] | alpha2, exp(-(lp2[n])/alpha2)) + 
                weibull_lpdf(YtYrdiff[n] | alpha3, exp(-(lp3[n])/alpha3));

    }
  }
}


