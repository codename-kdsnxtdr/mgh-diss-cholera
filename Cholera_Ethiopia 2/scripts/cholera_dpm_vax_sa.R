### cholera_dpm_vax_sa.R
#
# This script redefines aspects of the main transmission model script cholera_dpm_vax.R 
# to be used for sensitivity analysis. All main aspects of the model are the same, such as
# the odes and transmission parameters are the same. However, to improve run time and
# test vaccine related parameters, some embedded intervention functions were changed or
# removed:
#### hygiene_coverage function was removed
#### sanitation_coverage function was removed
#### cleanwater_coverage function was removed
#### forcing terms used to achieve a visual fit to the aseasonal peak were also removed 
#### vaccination function was changed so the coverage of a single or two-dose campaign could be varied 
#
# This script output the ode solutions over time given a single model run just as cholera_dpm_vax.R does.
# This function is sourced run_lhs_prcc.R to calculate outputs of interest for each iteration of sensitivity
# analysis, given LHS samples. 
#
###

### load relevant packages ----
library(deSolve)
library(readxl)


### define model function ----
cholera_dpm_vax_sa <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    
    # human population 
    P = S+V1+V2+E+Ev+I1+I2+I3+I4+Ia+Rc+Rcp+Ra 
    
    
    # WASH intervention 1: improve hygiene coverage 
    hyg_efct <- hyg_cov * 0.17/2
    HYG <- 1 - hyg_efct 
    
    # WASH intervention 2: improve sanitation coverage
    san_efct <- 7.7 * san_cov 
    SAN <- log10(K) - san_efct 
  
    # WASH intervention 3: improve clean water access
    cwa_efct <- cwa_cov * 1
    CWA <- 1 - cwa_efct
    
   
    # implement OCV campaigns ----
    vaccination <- function(t, v1_cov, rel_v2_cov) { # given a fixed OCV campaign adjust the coverage of one or two doses
      # define a fixed OCV campaign to begin at a random time
      vax_start = 2000 # arbitrary start time
      vax_dur1 = 7 # average from actual campaigns 
      vax_dur2 = 6 # average from actual campaigns
      vax_int = 72 # average from actual campaigns 
      
      # define time steps across whole simulation
      tps <- seq(0,365*10, by =1)
      
      # cumulative campaign rates holding vectors :: technically not necessary since fixed at one campaign, but fine
        # in case there are overlapping campaigns, because their rate would compound
      rho1_t <- numeric(length(tps)) # vaccination rate one dose
      rho2_t <- numeric(length(tps)) # vaccination rate two dose 
      
      # campaign specific rates holding vectors
      rho1_jt <- numeric(length(tps)) # vaccination rate one dose
      rho2_jt <- numeric(length(tps)) # vaccination rate two dose
      
      # for all timesteps where the campaign occurs set that value to one 
      rho1_jt[vax_start : (vax_start+vax_dur1)] <- 1
        
      # calculate the vaccination rate for one dose given the first-dose coverage and campaign duration
      rho1_j <- (-(log(1-(v1_cov))) / vax_dur1)
      # multiply vector of ones by the rate so that the rate is applied only at appropriate timesteps
      rho1_jt <- rho1_jt * rho1_j
        
      rho1_t <- rho1_t + rho1_jt # update the cumulative rate vector by adding the campaign-specific rate vector 
      
      # check if the v2 coverage is greater than zero, if not vaccination rate for two-doses will be zero
      if (rel_v2_cov > 0) {
        
        # same as above: for all timesteps where the campaign occurs set that value to one 
        rho2_jt[(vax_start+vax_dur1+vax_int): (vax_start+vax_dur1+vax_int+vax_dur2)] <- 1
        
        # make sure that V2 coverage is not greater than one :: technically not necessary since user inputs v2 coverage, but fine 
        cov2 <- rel_v2_cov
        cov2_capped <- ifelse(cov2 >= 1, 0.99, cov2)
        
        # calculate the vaccination rate for one dose given the first-dose coverage and campaign duration
        rho2_j <- (-(log(1-cov2_capped)) / vax_dur2)
        # multiply vector of ones by the rate so that the rate is applied only at appropriate timesteps
        rho2_jt <- rho2_jt * rho2_j
        
        rho2_t <- rho2_t + rho2_jt # update the cumulative rate vector by adding the campaign-specific rate vector 
            
        }
      
      # define timestep specific rates
      rho1 <- rho1_t[t]
      rho2 <- rho2_t[t]
      
      return(c("rho1" = rho1, "rho2" = rho2))
    }
    
    # run the function and extract values of interest
    vals <- vaccination(t, v1_cov, rel_v2_cov)
    rho1 <- vals["rho1"]
    rho2 <- vals["rho2"]
    
    
    # seasonality ----
    temp <- A_temp * sin((pi/6) * ((t/30+lag_tr) - phi_temp)) + temp_mean # temperature (degrees C) 
    rain <- A_rain * sin((pi/6) * ((t/30+lag_tr) - phi_rain)) + rain_mean # rainfall (mm)
    
    # calculated parameters ----
    beta_e <- beta_e0 * (1 + alpha_b * ((rain - rain_mean)/(rain_max - rain_mean)))
    
    KH <- KL/700 # IC50 of hyper-infectious bacteria, takes fewer hyper-infectious bacteria to make you sick 
    
    K_san <- 10^(SAN)
    
    lam_e <- beta_e * (BL/(KL + BL)) * CWA * HYG # environment-human FOI
    lam_h <- beta_h * (BH/(KH + BH)) * HYG # human-human FOI, impacted by hygiene intervention
    
    u <- zeta * (I1 + sigma_2*I2 + sigma_3*I3 + sigma_4*I4 + sigma_4*Ia) # human shedding rate of hyper-infectious bacteria into the environment
    n <- n_0 * (1 + alpha_n * ((temp - temp_mean)/(temp_max - temp_mean))) # intrinsic replication rate of environmental-state bacteria, impacted by sanitation intervention, impacted by water treatment intervention
    
    
    # odes ----
    dS <- mu*P + omega_a*Ra + omega_cp*Rcp + omega_v1*V1 + omega_v2*V2 - rho1*S - (lam_e+lam_h)*S - mu*S # susceptible
    
    dV1 <- rho1*S - omega_v1*V1 - rho2*V1 - (1-l1)^a*(lam_e+lam_h)*V1 - mu*V1 # vaccinated with one dose
    dV2 <- rho2*V1 - omega_v2*V2 - (1-l2)^a*(lam_e+lam_h)*V2 - mu*V2 # vaccinated with two doses
    
    dE <- (lam_e+lam_h)*S - (j+f+g)*gamma*E - (1-j-f-g)*gamma*E - mu*E # exposed
    
    dEv <- (1-l1)^a*(lam_e+lam_h)*V1 + (1-l2)^a*(lam_e+lam_h)*V2 - ((q*j)+(q*f)+(r*g))*gamma*Ev - (1-(q*j)-(q*f)-(r*g))*gamma*Ev - mu*Ev # exposed after vaccination
    
    dI1 <- j*gamma*(E + q*Ev) - (1-tau_t)*delta_1*theta_12*I1 - tau_t*theta_12*I1 - (1-tau_t)*(1-delta_1)*theta_12*I1 - mu*I1 # infected: severe case
    dI2 <- f*gamma*(E + q*Ev) + (1-tau_t)*(1-delta_1)*theta_12*I1 - (1-tau_t)*delta_2*theta_23*I2 - tau_t*theta_23*I2 - (1-tau_t)*(1-delta_2)*theta_23*I2 - mu*I2 # infected: moderate case
    dI3 <- g*gamma*(E + r*Ev) + (1-tau_t)*(1-delta_2)*theta_23*I2 - theta_34*I3 - mu*I3 # infected: mild case
    dI4 <- theta_34*I3 - tau_c*I4 - mu*I4 # infected: previously clinical asymptomatic case
    
    dIa <- (1-j-f-g)*gamma*E + (1-(q*j)-(q*f)-(r*g))*gamma*Ev - tau_a*Ia - mu*Ia # infected: asymptomatic case
    
    dRc <- tau_t*theta_12*I1 + tau_t*theta_23*I2 + tau_c*I4 - omega_cf*Rc - mu*Rc # recovered: clinical infection up to a year ago
    dRcp <- omega_cf*Rc - omega_cp*Rcp - mu*Rcp # recovered: clinical infection more than one year ago
    
    dRa <- tau_a*Ia - omega_a*Ra - mu*Ra # recovered: asymptomatic infection
    
    dBH <- u - b*v_H*BH - (1-b)*v_H*BH # hyperinfectious bacterial population
    dBL <- n*BL*(1-(BL/K)) + b*v_H*BH - v_L*BL # natural state bacterial population
    
    
    dCDc <- (1-tau_t)*delta_1*theta_12*I1 + (1-tau_t)*delta_2*theta_23*I2 # count compartment: cumulative cholera deaths
    dCInc <- (lam_e+lam_h)*S + (1-l1)^a*(lam_e+lam_h)*V1 + (1-l2)^a*(lam_e+lam_h)*V2 # count compartment: cumulative incidence
    
    # count compartment: cumulative reported deaths 
    dCRD <- ((rr_cd * (1-tau_t)*delta_1*theta_12*I1) +    # reported I1 cases --> cholera deaths
               (rr_cd * (1-tau_t)*delta_2*theta_23*I2))  # reported I2 cases --> cholera deaths
    
    # count compartment: cumulative reported cases 
    dCRInc <- ((rr_cd * (1-tau_t)*delta_1*theta_12*I1) +       # reported I1 cases --> cholera deaths 
                 (rr_t * tau_t*theta_12*I1) +                  # reported I1 cases --> treated 
                 (rr_cd * (1-tau_t)*delta_2*theta_23*I2) +     # reported I2 cases --> cholera deaths 
                 (rr_t * tau_t*theta_23*I2))                   # reported I2 cases --> treated
    
    # model outputs ----
    list(c(dS, dV1, dV2, dE, dEv, dI1, dI2, dI3, dI4, dIa, dRc, dRcp, dRa, dBH, dBL, dCDc, dCInc, dCRD, dCRInc))
  })
}


### initial conditions ----
initP <- 17700000
initV1 <- 0 
initV2 <- 0
initE <- 4011.49646455083 
initEv <- 0 
initI1 <- 112.331995953559
initI2 <- 288.210456951558 
initI3 <- 1045.0698452467
initI4 <- 340.1181209
initIa <- 1470.455851
initRc <- 1957401.12
initRcp <- 730176.988235497
initRa <- 268052.205848901
initS <- initP - initV1 - initV2 - initE - initEv - initI1 - initI2 - initI3 - initI4 - initIa - initRc - initRcp - initRa 

initBH <- 1.48e11
initBL <- 1.92e10

initCDc <- 0
initCInc <- 0
initCRD <- 0
initCRInc <- 0

istate <- c(S = initS,
            V1 = initV1, V2 = initV2,
            E = initE, Ev = initEv,
            I1 = initI1, I2 = initI2, I3 = initI3, I4 = initI4, Ia = initIa,
            Rc = initRc, Rcp = initRcp, Ra = initRa, 
            BH = initBH, BL = initBL,
            CDc = initCDc, CInc = initCInc,
            CRD = initCRD, CRInc = initCRInc)


### parameters ----
parameters <- list(
  mu = 1/(67.8*365), # human natural birth/death rate

  beta_e0 = 2.3e-4, # average environmental contact rate
  beta_h = 4.3e-6, # human contact rate

  gamma = 1/2, # incubation period
  j = 0.02, # proportion of cases severe, no vax
  f = 0.05, # proportion of cases moderate, no vax
  g = 0.18, #proportion of cases mild, no vax

  omega_cf = 1/3464.296, # first year loss of natural immunity, clinical case
  omega_cp = 1/1429.059, # after first year loss of natural immunity, clinical case
  omega_a = 1/(6*30), # loss of natural immunity, asymptomatic case

  theta_12 = 1/3, # untreated severe case improves to moderate
  theta_23 = 1/3, # untreated moderate case improves to mild
  theta_34 = 1/3, # untreated mild case improves to prev. clinical asymptomatic
  tau_c = 1, # prev. clinical asymptomatic case makes full recovery
  tau_a = 1, # non-clinical asymptomatic case makes full recovery

  tau_t = 0.8, # treatment rate of severe/moderate cases

  delta_1 = 0.5, # cholera death rate of untreated severe cases
  delta_2 = 0.1, # cholera death rate of untreated moderate cases

  zeta = 1e9, # severe case human bacterial shedding rate
  sigma_2 = 0.1, # relative shedding rate for moderate case
  sigma_3 = 0.01, # relative shedding rate for mild case
  sigma_4 = 0.0001, # relative shedding rate for asymptomatic case

  rr_cd = 0.29, # reporting rate of severe cholera cases who died from cholera
  rr_t = 1, # reporting rate of severe cholera cases who were treated

  v_H = 1, # hyperinfectious bacteria natural death rate/state-transition rate
  b = 0.1, # proportion of hyperinfectious bacteria secreted into the natural environment

  n_0 = 0.33, # mean natural state bacterial intrinsic proliferation
  v_L = 0.33, # natural bacteria natural death rate
  KL = 1e7, # concentration of V. cholerae in (natural environment) water that yields 50% chance of infection
  K = 1e10, # natural environment carrying capacity concentration

  temp_mean = 21.99, # mean temperature (celsius)
  temp_max = 23.98, # maximum temperature (celsius)
  rain_max = 1067245, # maximum rainfall (mm)
  A_temp = 1.18, # amplitude temperature equation
  A_rain = 319879.08, # amplitude rainfall equation
  phi_temp = 7.16, # phase shift temperature equation
  phi_rain = 8.87, # phase shift rainfall equation
  rain_mean = 414168.02, # mean rainfall (mm)

  alpha_n = 1.25, # natural state bacterial intrinsic proliferation sensitivity to temperature
  alpha_b = 2, # environmental contact rate sensitivity to rainfall
  lag_tr = 5, #temperature/rainfall lag

  omega_v1 = 1/(2*365), # 1-dose loss of vaccine-induced immunity
  omega_v2 = 1/(6.5*365), # 2-dose loss of vaccine-induced immunity
  a = 2, # vaccine dose assortitivity term
  l1 = 0.587, # 1-dose vaccine effectiveness
  l2 = 0.65, # 2-dose vaccine effectiveness
  q = 0.45, # vaccine risk reduction of severe+moderate cases
  r = 0.4, # vaccine risk reduction of mild cases

  # vaccination campaign coverage 
  v1_cov = 0, # one dose 
  rel_v2_cov = 0, # two dose (of those who got one dose, what proportion got a second)

  # baseline WASH coverage 
  hyg_cov = 0.05, # hygiene
  san_cov = 0.09, # sanitation
  cwa_cov = 0.5 # clean water
)


### timestep ----
time_start <- 1
time_stop <- 365 * 10
deltat <- 1 # 1 day

tps <- seq(time_start, time_stop, by = deltat)


### base model output (daily) ----
out_sa <- deSolve::ode(y = istate, times = tps, func = cholera_dpm_vax_sa, parms = parameters)
out_csa <- as.data.frame(out_sa)

out_csa
