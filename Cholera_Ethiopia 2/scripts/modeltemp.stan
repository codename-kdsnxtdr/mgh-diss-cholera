////// modeltemp.stan
//
// This script defines the information necessary for conducting Bayesian MCMC methods of 
// parameter estimation using rstan in the format that is recognizable by the rstan function:
// stan(). This information is necessary to estimate parameters values characteristic of a sine function
// representative of temperature data from Ethiopia between August 2022 and June 2024. 
// Key inputs include: the observed data, the model (sine function), priors and a likelihood funtion. 
//
// This script doesn't firectly generate outputs but is used by the stan() function
// in the 'run_modeltemp.R' script. 
//
//////


data { // inform Stan about what kind of data it will recieve
  int<lower=1> n_months; // a value representative of the number of data points
  vector[n_months] temps; // a vector 'temps' containing 'n_months' observations of national temperature data from Ethiopia
  vector[n_months] ts; // a time vector 'ts' containing as many timesteps as there are data points
}


parameters { // inform Stan about what parameters it will be sampling (parms of interest); provide boundaries, as needed
  real<lower=0> A_Te; // amplitude; should not be below zero
  real<lower=0, upper=12> phi_Te; // horizontal shift; bound within a single year to prevent confusion 
  real Te_mean; // mean temperature
  real sigma; // standard deviation 
}


transformed parameters { // inform Stan about how to transform parameter values to make it relevant to the data 
  vector[n_months] y; // define predicted temperature vector; should only be as long as there are data points
  y = (A_Te .* sin(pi()/6 .* (ts - phi_Te))) + Te_mean; // predict temperature given sampled parameters 
}


model {
  // define priors for paramters of interest
  A_Te ~ normal(2, 2); // mean chosen by looking at observed data
  phi_Te ~ normal(6, 3); // generic prior
  Te_mean ~ normal(22, 2); // mean chosen by looking at observed data
  sigma ~ uniform(-10, 10); // generic prior
  
  // define a likelihood function relating observed data with predicted temperature 
  temps ~ normal(y, sigma);
}


generated quantities {
  array[n_months] real pred_temp;
  // given predicted temperature distributions (y) for each time point,
  // reintroduce variability by selecting the reported temperature randomly from the posterior distribution 
  for (ii in 1:n_months) {
    pred_temp = normal_rng(y, sigma);
  }
}
