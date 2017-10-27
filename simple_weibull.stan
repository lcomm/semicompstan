data {
  // number of observations
  int<lower=0> N;
  // number of columns in design matrix, including intercept
  int<lower=1> P_1;
  // design matrix
  matrix[N, P_1] X_1;
  // observed event or censoring time
  real<lower=0> Y[N];
  // indicator of event observation
  int<lower=0,upper=1> dY[N];
}

parameters {
  // vector of regression parameters
  vector[P_1] beta1;
  
  // shape parameters (the one in exponent of time)
  // alpha > 1 -> hazard increases over time, more clumping
  real<lower=0> alpha1;
}

model {
  // linear predictors
  vector[N] lp1;
  lp1 = X_1 * beta1;
  
  // likelihood
  for (n in 1:N){
    if (dY[n] == 1) {
       Y[n] ~ weibull(alpha1, exp(-(lp1[n])/alpha1));
    } else {
      target += weibull_lccdf(Y[n] | alpha1, exp(-(lp1[n])/alpha1));
    }
  }
  
}
