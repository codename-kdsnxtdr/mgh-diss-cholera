////// modelrain.stan
//
// This script defines the infomration necessary for conducting Bayesian MCMC methods of 
// parameter estimation using rstan in the format that is recognizable by the rstan function:
// stan(). This information is necessary to estimate parameters values characteristic of a sine function
// representative of rainfall data from Ethiopia between August 2022 and June 2024. 
// Key inputs include: the observed data, the model (sine function), priors and a likelihood funtion. 
//
// This script doesn't firectly generate outputs but is used by the stan() function
// in the 'run_modelrain.R' script. 
//
//////


data { // inform Stan about what kind of data it will recieve
  int<lower=1> n_months; // a value representative of the number of data points
  vector[n_months] rain; // a vector 'temps' containing 'n_months' observations of national rainfall data from Ethiopia
  vector[n_months] ts; // a time vector 'ts' containing as many timesteps as there are data points
}


parameters { // inform Stan about what parameters it will be sampling (parms of interest); provide boundaries, as needed
  real<lower=0> A_RaC; // amplitude/100000; should not be below zero
  real<lower=0, upper=12> phi_Ra;  // horizontal shift; bound within a single year to prevent confusion 
  real Ra_meanC; // mean rainfall divided by 100000
  real sigma; // standard deviation 
}


transformed parameters { // inform Stan about how to transform parameter values to make it relevant to the data 
  vector[n_months] y; // define predicted rainfall vector; should only be as long as there are data points
  y = A_RaC .* sin(pi()/6 .* (ts - phi_Ra)) + Ra_meanC; // predict rainfall/100000 given sampled parameters 
  
  // actual values of rainfall in mm was divided by 100,000 to prevent confusion/complexity for the algorithm and make definition of priors simpler
  // values were converted back to mm here to be matched with the data 
  real A_Ra = A_RaC * 100000; // convert back to mm
  real Ra_mean = Ra_meanC * 100000; // convert back to mm
}


model {
  // // define priors for paramters of interest
  A_RaC ~ normal(5, 2); // mean chosen by looking at observed data
  phi_Ra ~ normal(6, 3); // generic prior
  Ra_meanC ~ normal(5, 2); // mean chosen by looking at observed data
  sigma ~ uniform(-10, 10); // generic prior
  
  // define a likelihood function relating observed data with predicted rainfall 
  rain ~ normal(y, sigma);
}


generated quantities {
  array[n_months] real pred_rain;
  // given predicted temperature distributions (y) for each time point, covert to mm and
  // reintroduce variability by selecting the reported rainfall randomly from the posterior distribution 
  for (ii in 1:n_months) {
    pred_rain = normal_rng(y * 100000, sigma);
  }
}
