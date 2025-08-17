### run_lhs_prcc.R
#
# This script is used to produce outputs from sensitivity analysis using the PRCC
# method applied to parameters sampled generated via LHS.
#
# The script defines a function that calculates model outputs of interest using the 
# 'cholera_dpm_vax_sa.R' script. It also produces a matrix of LHS samples from provided
# parameter distributions. These samples, in a for loop are fed to the function 
# generating a matrix containing all LHS samples with thir associated model outputs. 
# Finally this matrix is fed to the epi.prcc() funtion from the epiR package to calculate 
# partial rank correlation coefficients.
#
###

### load relevant packages ----
library(sensitivity)
library(epiR)
library(dplyr)
library(lhs)
library(readxl)
library(deSolve)


### source modified trnamsission model that will be used to calculate key outputs
source("cholera_dpm_vax_sa.R")


### define a function that will take in a list of parameter values and return key outputs of the transmission model
#
cholera_dpm_vax_safun <- function(parameters) {
  # produce out object containing transmission model solutions across all timesteps given defined parameters
  with(as.list(c(parameters)), {
    out <- tryCatch({ 
      # provide instructions for what to do if a specific parameters sets results in warnings or solution errors
      withCallingHandlers(
        {ode(y = istate, times = tps, func = cholera_dpm_vax_sa, parms = parameters)},
      warning = function(w) { # if there is a warning, keep going without printing
        invokeRestart("muffleWarning")}
      )}, 
      error = function(e) { # if there is an error return out as NULL
        return(NULL)}
    )
    
    # if out indicates a solving error (i.e. is NULL), return outputs of interest as NA
    if (is.null(out)) {
      return(c(NA_real_, NA_real_))
    }
    
    out_c <- as.data.frame(out) 
    
    # if the solver wasn't able to execute all timesteps, return outputs of interest as NA
    if (!all(tps %in% out_c$time)) {
      return(c(NA_real_, NA_real_))
    }
    
    # if any value of model outputs required to calculate outputs of interest are infinite or NA, return outputs of interest as NA
    if (any(is.na(out_c$CInc)) || 
        any(is.infinite(out_c$CInc)) ||
        any(is.na(out_c$CDc)) ||
        any(is.infinite(out_c$CInc))) {
      return(c(NA_real_, NA_real_))
    }
    
    # if the run is successful, calculate outputs of interest
    c_inc_2528 <- out_c$CInc[3270] - out_c$CInc[2190] # total case incidence between Jan2025 and Jan2028
    d_inc_2528 <- out_c$CDc[3270] - out_c$CDc[2190] # deaths between Jan2025 and Jan2028
    
    cr_inc_2528 <- out_c$CRInc[3270] - out_c$CRInc[2190] # for calc CFR: reported case incidence between Jan2025 and Jan2028 
    dr_inc_2528 <- out_c$CRD[3270] - out_c$CRD[2190] # for calc CFR: reported deaths between Jan2025 and Jan2028  
    cfr_2528 <- dr_inc_2528 / cr_inc_2528 # CFR between Jan2025 and Jan2028
    
    return(c(c_inc_2528,d_inc_2528,cfr_2528)) # return vector of outputs of interest
  })
}


### define parameter ranges that LHS will sample from
parm.ranges = list(
  # varied parameters for sensitivity analysis
  beta_e0 = c(1.15e-4,4.6e-4), # average environmental infectious contact rate
  beta_h = c(2.15e-6,8.6e-6), # human infectious contact rate
  gamma = c(1/5, 1/0.5), # incubation period 
  omega_cp = c(1/730,1/1460), # after first year loss of natural immunity, clinical case 
  omega_a = c(1/(9*30), 1/(6*30)), # loss of natural immunity, asymptomatic case
  tau_a = c(1/2,1), # non-clinical asymptomatic case makes full recovery
  tau_t = c(0.73,0.88), # recovery of treated severe cases
  delta_2 = c(0.1,0.25), # cholera death rate of untreated moderate cases
  zeta = c(1e9,1e12), # severe case human bacterial shedding rate
  K = c(1e2,1e10), # natural environment carrying capacity concentration
  alpha_n = c(1,1.8), # natural state bacterial intrinsic proliferation sensitivity to temperature
  alpha_b = c(1,2), # environmental contact rate sensitivity to rainfall
  omega_v1 = c(1/(2*365),1/(1*365)), # loss of vaccine-induced immunity 1D
  omega_v2 = c(1/(7*365),1/(5*365)), # loss of vaccine-induced immunity 2D
  a = c(1,3), # vaccine assortitivity term
  
  # not varied parameters
  mu = c(1/(67.8*365),1/(67.8*365)), # human natural birth/death rate
  j = c(0.02,0.02), # proportion of cases severe, no vax
  f = c(0.05,0.05), # proportion of cases moderate, no vax
  g = c(0.18,0.18), # proportion of cases mild, no vax
  omega_cf = c(1/3464.296,1/3464.296), # first year loss of natural immunity, clinical case
  theta_12 = c(1/3,1/3), # untreated severe case improves to moderate
  theta_23 = c(1/3,1/3), # untreated moderate case improves to mild
  theta_34 = c(1/3,1/3), # untreated mild case improves to prev. clinical asymptomatic 
  tau_c = c(1,1), # prev. clinical asymptomatic case makes full recovery
  delta_1 = c(0.5,0.5), # cholera death rate of untreated severe case
  sigma_2 = c(0.1,0.1), # relative shedding rate for moderate case
  sigma_3 = c(0.01,0.01), # relative shedding rate for mild case 
  sigma_4 = c(0.0001,0.0001), # relative shedding rate for asymptomatic case
  rr_cd = c(0.29,0.29), # reporting rate of severe cholera cases who died from cholera
  rr_t = c(1,1), # reporting rate of severe cholera cases who were treated
  v_H = c(5,5), # hyperinfectious bacteria natural death rate/state-transition rate
  b = c(0.1,0.1), # proportion of hyperinfectious bacteria secreted into the natural environment
  n_0 = c(0.33,0.33), # mean natural state bacterial intrinsic proliferation 
  v_L = c(0.33,0.33), # natural bacteria natural death rate
  KL = c(1e7,1e7), # concentration of hyperinfectious V. cholerae in (natural environment) water that yields 50% chance of infection
  temp_mean= c(21.99,21.99), # mean temperature (celsius)
  temp_max = c(23.98,23.98), # maximum temperature (celsius)
  rain_max = c(1067245,1067245), # maximum rainfall (mm)
  A_temp = c(1.18,1.18), # amplitude temperature equation
  A_rain = c(319879.08,319879.08), # amplitude rainfall equation
  phi_temp = c(7.16,7.16), # phase shift temperature equation
  phi_rain = c(8.87,8.87), # phase shift rainfall equation
  rain_mean = c(414168.02,414168.02), # mean rainfall (mm)
  lag_tr = c(5,5), # temperature and rainfall lag
  l1 = c(0.587,0.587), # vaccine efficacy 1D
  l2 = c(0.65,0.65), # vaccine efficacy 2D
  q = c(0.45,0.45), # vaccine risk reduction of severe+moderate cases
  r = c(0.4,0.4),  # vaccine risk reduction of mild cases
  v1_cov = c(0,0), # vaccination campaign coverage: one-dose 
  rel_v2_cov = c(0,0), # vaccination campaign coverage: two-doses
  hyg_cov = c(0.05,0.05), # baseline WASH coverage: hygiene
  san_cov = c(0.09,0.09), # baseline WASH coverage: sanitation
  cwa_cov = c(0.5,0.5) # baseline WASH coverage: clean water access
)


### produce sensitivity analysis outputs 
num_parms <- 53 # 15 of which are  being varied 
num_samples <- 1500 # 100 samples/varied parameter

lhs_matrix <- randomLHS(num_samples, num_parms) # construct LHS matrix (not parameter specific yet)

parm.samples <- matrix(NA, nrow = num_samples, ncol = num_parms) # construct an empty matrix to house all samples for all paramters
colnames(parm.samples) <- names(parm.ranges) # give the columns parameter names

# apply the multiplier of the constructed random LHS matrix to fit the specific parameter range 
for (i in seq_len(num_parms)) { 
  min_val <- parm.ranges[[i]][1]
  max_val <- parm.ranges[[i]][2]
  # increase the minimum parameter value by some proportion within the defined range, according to LHS matrix
  parm.samples[,i] <- lhs_matrix[, i] * (max_val - min_val) + min_val 
}


parm.samples.df <- as.data.frame(parm.samples)
# add columns for the outputs of interest, to store the results for every sample
parm.samples.df[["c_inc_2528"]] <- 0.0
parm.samples.df[["d_inc_2528"]] <- 0.0
parm.samples.df[["cfr_2528"]] <- 0.0



f_name <- "lhs.prcc.sims.rds" # define a save file name
if (!file.exists(f_name)) { # prevent overwriting a prexisting save file
  start_time <- Sys.time()
  # for every sample, use above function to calculate outputs of interest and then store the values in the dataframe from above
  for (rr in 1:nrow(parm.samples.df)) {
    pp <- parm.samples.df[rr,]
    parm.samples.df[rr, "c_inc_2528"] <- cholera_dpm_vax_safun(pp)[1]
    parm.samples.df[rr, "d_inc_2528"] <- cholera_dpm_vax_safun(pp)[2]
    parm.samples.df[rr, "cfr_2528"] <- cholera_dpm_vax_safun(pp)[3]
  }
  end_time <- Sys.time()
  saveRDS(parm.samples.df, file = f_name) # save data frame as an external file
}

lhs.sims <- readRDS(f_name) # load in saved file



### calculate PRCC
parm.names <- setdiff(names(lhs.sims), c("c_inc_2528", "d_inc_2528", "cfr_2528")) # differentiate outputs from samples

# remove all columns for parameters that weren't varied; isoloate varied parameters
for (p in parm.names) {
  if (length(unique(lhs.sims[[p]])) == 1) {
    lhs.sims <- lhs.sims %>% select(-all_of(p))
  }
}

varied_parm.names <- setdiff(names(lhs.sims), c("c_inc_2528", "d_inc_2528", "cfr_2528")) # differentiate outputs from samples, varied only

# remove all failed runs, where outputs were returned as NA and produce output specific data frames
prcc_df_cinc <- lhs.sims[!is.na(lhs.sims$c_inc_2528), c(varied_parm.names, "c_inc_2528")] 
prcc_df_dinc <- lhs.sims[!is.na(lhs.sims$c_inc_2528), c(varied_parm.names, "d_inc_2528")]

# calculate separate prccs for each output of interest
prcc_result_cinc <- epi.prcc(prcc_df_cinc, sided.test = 2)
prcc_result_dinc <- epi.prcc(prcc_df_dinc, sided.test = 2)


