### fit_cholera_dpm_vax.R
#
# This script takes in relevant data and formats it for use by modeltemp.stan to
# execute parameter estimation using Bayesian MCMC methods. Parameters being estimated
# were used in a sine function representative of temperature over time fitted to national
# Ethiopian temperature data.
#
# This script outputs the numerical posterior outputs of the Stan model for parameters of
# interest and plots that can be used to validate the sufficiency of these estimates including, 
# posterior distributions across chains, trace plots and model outputs plotted against observed data.
#
###


### load relevant packages----
library(rstan)
library(bayesplot)
library(ggplot2)
library(esquisse)
library(tibble)
library(tidyr)
library(gridExtra)
library(readxl)


### run the Stan model ----
## set Stan options
rstan_options(auto_write = TRUE) # prevent recompilation to C++ to save time in subsequent runs
options(mc.cores = parallel::detectCores()) # enable parallel execution of  MCMC chains to reduce computation time

## format raw data for recognition by Stan model
dat_raw <- read_excel("data/raw_data/Cholera_Data_Monthly_Mule_1.xlsx") # load cholera incidence data
rinc_dat <- dat_raw$incidence[6:length(dat$incidence)] # exclude first 5 points

## define fixed parameters
n_days_sim <- 10*365 # number of simulation days 
n_months_sim <- floor(n_days_sim/30)
burnin_months <- 4*12 # months to run before fitting to data 
n_months_data <- length(rinc_dat) # months of data


mu  <- 1/(67.8*365) 
gamma <- 1/2 # incubation period 
j <- 0.02 # proportion of cases severe, no vax
f <- 0.05 # proportion of cases moderate, no vax
g <- 0.18 #proportion of cases mild, no vax 
omega_cf <- 1/3464.296 # first year loss of natural immunity, clinical case
omega_cp <- 1/1429.059 # after first year loss of natural immunity, clinical case 
omega_a <- 1/(6*30) # loss of natural immunity, asymptomatic case
theta_12 <- 1/3 # untreated severe case improves to moderate
theta_23 <- 1/3 # untreated moderate case improves to mild
theta_34 <- 1/3 # untreated mild case improves to prev. clinical asymptomatic 
tau_c <- 1 # prev. clinical asymptomatic case makes full recovery
tau_a <- 1 # non-clinical asymptomatic case makes full recovery
tau_t <- 0.8 #0.85875, #0.613, # treatment rate of severe/moderate cases
delta_1 <- 0.5 # cholera death rate of untreated severe cases
delta_2 <- 0.1#0.25, # cholera death rate of untreated moderate cases
zeta <- 1e9 # severe case human bacterial shedding rate
sigma_2 <- 0.1 # relative shedding rate for moderate case
sigma_3 <- 0.01 # relative shedding rate for mild case 
sigma_4 <- 0.0001 # relative shedding rate for asymptomatic case
rr_cd <- 0.29 #0.43, # reporting rate of severe cholera cases who died from cholera
rr_t <- 1 # reporting rate of severe cholera cases who were treated
v_H <- 1 # hyperinfectious bacteria natural death rate/state-transition rate
b <- 0.1 # proportion of hyperinfectious bacteria secreted into the natural environment
n_0 <- 0.33 # mean natural state bacterial intrinsic proliferation 
v_L <- 0.33 # natural bacteria natural death rate
KL <- 1e7 # concentration of V. cholerae in (natural environment) water that yields 50% chance of infection
K <- 1e10 # natural environment carrying capacity concentration
temp_mean <- 21.99 #21.97583, # mean temperature (celsius)
temp_max <- 23.98 # maximum temperature (celsius)
rain_max <- 1067245# maximum rainfall (mm)
A_temp <- 1.18 #1.11, # amplitude temperature equation
A_rain <- 319879.08 #337062.96, # amplitude rainfall equation
phi_temp <- 7.16 #2.21, # phase shift temperature equation
phi_rain <- 8.87 # phase shift rainfall equation
rain_mean <- 414168.02# mean rainfall (mm)
lag_tr <- 5#temperature/rainfall lag
omega_v1  <- 1/(2*365) # 1-dose loss of vaccine-induced immunity
omega_v2 <- 1/(5*365) # 2-dose loss of vaccine-induced immunity
a <- 2 # vaccine dose assortativity term
l1 <- 0.587 # 1-dose vaccine effectiveness 
l2 <- 0.65 # 2-dose vaccine effectiveness
q <- 0.45 # vaccine risk reduction of severe+moderate cases
r <- 0.4 # vaccine risk reduction of mild cases
hyg_0 <- 0.05
san_0 <- 0.07
cwa_0 <- 0.50
rho1 <- 0
rho2 <- 0



## define initial conditions
P0 <- 17700000
y0 <- c(
  S = (P0 - 4011.49646455083 - 112.331995953559 - 288.210456951558 - 1045.0698452467 - 340.1181209 - 1470.455851 - 1957401.12 - 730176.988235497 - 268052.205848901),
  V1 = 0,
  V2 = 0,
  E = 4011.49646455083,
  Ev = 0,
  I1 = 112.331995953559,
  I2 = 288.210456951558,
  I3 = 1045.0698452467,
  I4 = 340.1181209,
  Ia = 1470.455851,
  Rc = 1957401.12,
  Rcp = 730176.988235497,
  Ra = 268052.205848901
)
y0["BH"] <- 0
y0["BL"] <- 0
y0["CRInc"] <- 0


## set initial coditions for sampling guide convergence 
init_pars <-function() { # initial conditions for parameters selected from a narrow distribution
  list(
    beta_e0 = 2.15e-4,
    beta_h = 4.4e-6,
    alpha_b = 2,
    alpha_n = 1.25
  )
}


## compile parameters for Stan data block
data <- list(
  n_days_sim = n_days_sim,
  n_months_sim = n_months_sim,
  n_months_data = n_months_data,
  burnin_months = burnin_months, 
  mu = mu,       
  gamma = gamma, 
  j = j, 
  f = f,
  g = g, 
  omega_cf = omega_cf, 
  omega_cp = omega_cp,
  omega_a = omega_a,
  theta_12 = theta_12, 
  theta_23 = theta_23,  
  theta_34 = theta_34,
  tau_c = tau_c,
  tau_a = tau_a, 
  tau_t = tau_t, 
  delta_1 = delta_1, 
  delta_2 = delta_2, 
  zeta = zeta, 
  sigma_2 = sigma_2, 
  sigma_3 = sigma_3,
  sigma_4 = sigma_4, 
  rr_cd = rr_cd, 
  rr_t = rr_t, 
  v_H = v_H,
  b = b, 
  n_0 = n_0,
  v_L = v_L,
  KL = KL,
  K = K,
  temp_mean = temp_mean, 
  temp_max = temp_max, 
  rain_max = rain_max, 
  A_temp = A_temp, 
  A_rain = A_rain, 
  phi_temp = phi_temp,
  phi_rain = phi_rain,
  rain_mean = rain_mean, 
  lag_tr = lag_tr,
  omega_v1 = omega_v1,
  omega_v2 = omega_v2,
  a = a, 
  l1 = l1, 
  l2 = l2, 
  q = q,
  r = r, 
  hyg_0 = hyg_0,
  san_0 = san_0,
  cwa_0 = cwa_0, 
  rho1 = rho1,
  rho2 = rho2, 
  y0 = y0,
  rinc_dat = rinc_dat,
  compute_likelihood = 1
)


## define MCMC iterations and number of chains
niter <- 1001 #4000
nwarmup <- 1000 #2000
nchains <- 1 #4


## run Stan model
cholera_dpm_vax_samples <- stan(
  "scripts/fit_cholera_dpm_vax.stan", # call Stan script
  iter = niter,
  warmup = nwarmup,
  data = data,
  chains = nchains,
  cores = nchains,
  seed = 31415, # set seed for reproducability
  refresh = 100,
  init = replicate(nchains, init_pars(), simplify = FALSE)
)



### produce outputs from Stan model run [NO OUTPUTS PRODUCED] ----

pars <- c('beta_e0', 'beta_h', 'alpha_n', 'alpha_b') # define modeled/transformed parameters of interest
print(cholera_dpm_vax_samples, pars = pars) # print model output :: estimations for parameter values
