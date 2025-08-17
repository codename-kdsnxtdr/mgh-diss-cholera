### cholera_dpm_vax.R
#
# This script takes in timestep, parameters, initial condition and a cholera model function to 
# solve the system of odes contained by the model function, representing cholera transmission.
#
# This script outputs a data frame of compartmental values of the model at each daily timestep and a 
# data frame containing the observed and the calculated monthly reported case incidence, deaths and cfr
# calculated by an external function, given the daily moddel output
#
###


### load relevant packages ----
library(deSolve) # to solve model system of odes
library(readxl) # to read in data when calculating model outputs


### define model function ----
cholera_dpm_vax <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    
    # human population 
    P = S+V1+V2+E+Ev+I1+I2+I3+I4+Ia+Rc+Rcp+Ra 
    
    # WASH intervention 1: improve hygiene coverage ----
    # simulate hygiene coverage improvement campaign: given a starting coverage, this function improves coverage between specified start and end points, assuming linearity
    hygiene_coverage <- function(t, hyg_starts, hyg_ends, hyg_goals, hyg_0) {
      current <- hyg_0
      if (length(hyg_starts) == 0) { # if no campaigns specified, return starting coverage
        current <- hyg_0
      } else {
        for (j in seq_along(hyg_starts)) { # carry out improvements for each index of the start vector (per campaign)
          result <- numeric(length(hyg_ends[j] - hyg_starts[j])) # results vector is the length of the campaign
          result[1] <- current # first index of the result = starting coverage
          if (t < hyg_starts[j]) { # if current timestep before campaign period, no change to coverage
            current <- current
          } else if (t >= hyg_starts[j] && t <= hyg_ends[j]) { # if current timestep in the campaign period, extrapolate linearly to the campaign end
            x <- t - hyg_starts[j] # keep track of distance from campaign start for slope calculation
            for (i in seq_along(x)) {
              m <- (hyg_goals[j] - result[1]) / (hyg_ends[j] - hyg_starts[j]) # calculate linear slope between start cov and goal coverage over desired timeframe
              result[i] <- pmin(1, current + m*x[i]) # calculate new coverage based on distance from start time
              
              current <- result[i]
            }
          } else { # if current timestep after campaign ends, maintain the goal coverage 
            current <- hyg_goals[j]
          }
        }
      }
      return(current)
    }
    
    hval <- hygiene_coverage(t, hyg_starts, hyg_ends, hyg_goals, hyg_0) # run the function
    hyg_cov <- hval[1] # extract function output
    
    hyg_efct <- hyg_cov * 0.17/2 # calculate effectiveness, assuming some combination of effectiveness and adherence 
    
    # aseasonal outbreak peak calibration: HYG forcing term 
    if (t >= 1930 && t < 2220) { # month 64 to 74
      hyg_efct <- 0.40
    }
    
    HYG <- 1 - hyg_efct # calculate transmission parameter reductive factor 
   
    
    # WASH intervention 2: improve sanitation coverage ----
    # simulate sanitation coverage improvement campaign: given a starting coverage, this function improves coverage between specified start and end points, assuming linearity
    sanitation_coverage <- function(t, san_starts, san_ends, san_goals, san_0) {
      current <- san_0
      if (length(san_starts) == 0) { # if no campaigns specified, return starting coverage
        current <- san_0
      } else {
        for (j in seq_along(san_starts)) { # carry out improvements for each index of the start vector (per campaign)
          result <- numeric(length(san_ends[j] - san_starts[j])) # results vector is the length of the campaign
          result[1] <- current  # first index of the result = starting coverage
          if (t < san_starts[j]) { # if current timestep before campaign period, no change to coverage
            current <- current
          } else if (t >= san_starts[j] && t <= san_ends[j]) { # if current timestep in the campaign period, extrapolate linearly to the campaign end
            x <- t - san_starts[j] # keep track of distance from campaign start for slope calculation
            for (i in seq_along(x)) {
              m <- (san_goals[j] - result[1]) / (san_ends[j] - san_starts[j]) # calculate linear slope between start cov and goal coverage over desired timeframe
              result[i] <- pmin(1, current + m*x[i]) # calculate new coverage based on distance from start time
              
              current <- result[i]
            }
          } else { # if current timestep after campaign ends, maintain the goal coverage 
            current <- san_goals[j]
          }
        }
      }
      return(current)
    }
    
    sval <- sanitation_coverage(t, san_starts, san_ends, san_goals, san_0) # run the function
    san_cov <- sval[1] # extract function output
    
    # aseasonal outbreak peak calibration: SAN forcing term 
    if (t >= 1930 && t < 2220) { # month 64 to 74
      san_cov <- 0.4
    }
    
    # define intervention effect on carrying capacity K
    san_efct <- 7.7 * san_cov # coverage determine the log-scale reduction of K
    SAN <- log10(K) - san_efct 
    
    
    # WASH intervention 3: improve clean water access ----
    # simulate clean water access improvement campaign: given a starting coverage, this function improves coverage between specified start and end points, assuming linearity
    cleanwater_coverage <- function(t, cwa_starts, cwa_ends, cwa_goals, cwa_0) {
      current <- cwa_0
      if (length(cwa_starts) == 0) { # if no campaigns specified, return starting coverage
        current <- cwa_0
      } else {
        for (j in seq_along(cwa_starts)) { # carry out improvements for each index of the start vector (per campaign)
          result <- numeric(length(cwa_ends[j] - cwa_starts[j])) # results vector is the length of the campaign
          result[1] <- current  # first index of the result = starting coverage
          if (t < cwa_starts[j]) { # if current timestep before campaign period, no change to coverage
            current <- current
          } else if (t >= cwa_starts[j] && t <= cwa_ends[j]) { # if current timestep in the campaign period, extrapolate linearly to the campaign end
            x <- t - cwa_starts[j] # keep track of distance from campaign start for slope calculation
            for (i in seq_along(x)) {
              m <- (cwa_goals[j] - result[1]) / (cwa_ends[j] - cwa_starts[j]) # calculate linear slope between start cov and goal coverage over desired timeframe
              result[i] <- pmin(1, current + m*x[i]) # calculate new coverage based on distance from start time
              
              current <- result[i]
            }
          } else { # if current timestep after campaign ends, maintain the goal coverage 
            current <- cwa_goals[j]
          }
        }
      }
      return(current)
    }
    
    cval <- cleanwater_coverage(t, cwa_starts, cwa_ends, cwa_goals, cwa_0) #run the function
    cwa_cov <- cval[1] # extract function output
    
    cwa_efct <- cwa_cov * 1
    
    # aseasonal calibration: CWA forcing term
    if (t >= 1930 && t < 2220) { # month 64 to 74
      cwa_efct <- 0.7
    }
    
    CWA <- 1 - cwa_efct # calculate transmission reductive factor 
    
    
    # implement OCV campaigns ----
    # provided vectors containing index-wise OCV campaign information, produce vectors containing the aggregate vaccination rate for one- and two-doses at each timepoint
    vaccination <- function(t, vax_starts, vax2D, vax_dur1, pvax1, vax_dur2, vax_int, pvax2) {
      tps <- seq(0,365*10, by =1) # define timepoints across full simulation 
      rho1_t <- numeric(length(tps)) # initialize final aggregate vector for rho1
      rho2_t <- numeric(length(tps)) # initialize final aggregate vector for rho2
      
      if (length(vax_starts) == 0 ) { # if vax_starts is empty, no campaings so rates are zero
        return(c("rho1" = 0, "rho2" = 0))
      } else { # if campaigns...
        for (j in seq_along(vax_starts)) { # loop through each campaign and calculate campaign rho1 and rho2 vevtors 
          rho1_jt <- numeric(length(tps)) # initialize campaign-specific vector for rho1
          rho2_jt <- numeric(length(tps)) # initialize campaign-specific vector for rho2
          
          rho1_jt[vax_starts[j] : (vax_starts[j]+vax_dur1[j])] <- 1 # for all timesteps rho1 should be on, vector index becomes one, else 0
          
          rho1_j <- (-(log(1-(pvax1[j]/initP))) / vax_dur1[j]) # calculate coverage and duration-dependent rate for rho1
          rho1_jt <- rho1_jt * rho1_j # multiply calculated rate to vector of 0s and 1s, so at all times rho1 should be on value is equal to calulated rho1
          
          rho1_t <- rho1_t + rho1_jt # add campaign specific vector to aggregate vector
          
          if (vax2D[j] == TRUE) { # if campaign is two-dose format, calculate rho2, else rho2 equals 0
            
            rho2_jt[(vax_starts[j]+vax_dur1[j]+vax_int[j]): (vax_starts[j]+vax_dur1[j]+vax_int[j]+vax_dur2[j])] <- 1 # for all timesteps rho2 should be on, vector index becomes one, else 0
            
            cov2 <- pvax2[j]/V1 # calculate coverage relative to those who received the first dose
            cov2_capped <- ifelse(cov2 >= 1, 0.99, cov2) # prevent log(0) rate calculation failure
            
            rho2_j <- (-(log(1-cov2_capped)) / vax_dur2[j]) # calculate coverage and duration-dependent rate for rho2
            rho2_jt <- rho2_jt * rho2_j # multiply calculated rate to vector of 0s and 1s, so at all times rho1 should be on value is equal to calulated rho2
            
            rho2_t <- rho2_t + rho2_jt # add campaign specific vector to aggregate vector
            
          }
        }
      }
      # extract aggregate rho values for the current timestep
      rho1 <- rho1_t[t]
      rho2 <- rho2_t[t]
      
      return(c("rho1" = rho1, "rho2" = rho2))
    }
    
    vals <- vaccination(t, vax_starts, vax2D, vax_dur1, pvax1, vax_dur2, vax_int, pvax2) # run the function
    
    rho1 <- vals["rho1"] # extract function outputs
    rho2 <- vals["rho2"] # extract function outputs
    
    
    # seasonality ----
    temp <- A_temp * sin((pi/6) * ((t/30+lag_tr) - phi_temp)) + temp_mean # temperature (degrees C) 
    rain <- A_rain * sin((pi/6) * ((t/30+lag_tr) - phi_rain)) + rain_mean # rainfall (mm)
    
    # calculated parameters ----
    beta_e <- beta_e0 * (1 + alpha_b * ((rain - rain_mean)/(rain_max - rain_mean))) # rainfall-dependent envrionmental infectious contact rate
    
    # aseasonal calibration: infectious contact forcing terms
    if (t >= 1790 && t < 1930) { # month 60 to 64.6
      beta_e <- beta_e + beta_e0 * 1.7
      beta_h <- beta_h * 5
    }
    
    KH <- KL/700 # hyper-infectious bacteria ID50 relative to natural-state bacteria ID50
    K_san <- 10^(SAN) # environmental carrying capacity of natural-state bacteria given sanitation intervention effect
    
    lam_e <- beta_e * (BL/(KL + BL)) * CWA * HYG # environment-human FOI, given hygiene and clean water access interventional effects 
    lam_h <- beta_h * (BH/(KH + BH)) * HYG # human-human FOI, give hygiene intervention effect
    
    u <- zeta * (I1 + sigma_2*I2 + sigma_3*I3 + sigma_4*I4 + sigma_4*Ia) # human shedding rate of hyper-infectious bacteria into the environment
    n <- n_0 * (1 + alpha_n * ((temp - temp_mean)/(temp_max - temp_mean))) # intrinsic replication rate of natural-state bacteria
    
    # CFR calibration: treatment rate forcing term 
    if (t >= 1828 && t <= 3650) { # month 61 (jan2024) to 64.6 
      tau_t <- 0.90
    }
    
    
    # ordinary differential equations ----
    dS <- mu*P + omega_a*Ra + omega_cp*Rcp + omega_v1*V1 + omega_v2*V2 - rho1*S - (lam_e+lam_h)*S - mu*S # susceptible
    
    dV1 <- rho1*S - omega_v1*V1 - rho2*V1 - (1-l1)^a*(lam_e+lam_h)*V1 - mu*V1 # vaccinated with one-dose OCV
    dV2 <- rho2*V1 - omega_v2*V2 - (1-l2)^a*(lam_e+lam_h)*V2 - mu*V2 # vaccinated with two-dose OCV
    
    dE <- (lam_e+lam_h)*S - (j+f+g)*gamma*E - (1-j-f-g)*gamma*E - mu*E # exposed
    
    dEv <- (1-l1)^a*(lam_e+lam_h)*V1 + (1-l2)^a*(lam_e+lam_h)*V2 - ((q*j)+(q*f)+(r*g))*gamma*Ev - (1-(q*j)-(q*f)-(r*g))*gamma*Ev - mu*Ev # exposed with previous vaccination
    
    dI1 <- j*gamma*(E + q*Ev) - (1-tau_t)*delta_1*theta_12*I1 - tau_t*theta_12*I1 - (1-tau_t)*(1-delta_1)*theta_12*I1 - mu*I1 # severe infection
    dI2 <- f*gamma*(E + q*Ev) + (1-tau_t)*(1-delta_1)*theta_12*I1 - (1-tau_t)*delta_2*theta_23*I2 - tau_t*theta_23*I2 - (1-tau_t)*(1-delta_2)*theta_23*I2 - mu*I2 # moderate infection
    dI3 <- g*gamma*(E + r*Ev) + (1-tau_t)*(1-delta_2)*theta_23*I2 - theta_34*I3 - mu*I3 # mild infection
    dI4 <- theta_34*I3 - tau_c*I4 - mu*I4 # previously clinical asymptomatic infection
    
    dIa <- (1-j-f-g)*gamma*E + (1-(q*j)-(q*f)-(r*g))*gamma*Ev - tau_a*Ia - mu*Ia # asymptomatic infection
    
    dRc <- tau_t*theta_12*I1 + tau_t*theta_23*I2 + tau_c*I4 - omega_cf*Rc - mu*Rc # recovered from clinical infection, within a year-ish
    dRcp <- omega_cf*Rc - omega_cp*Rcp - mu*Rcp # recovered from clinical infection, after a year-ish
    
    dRa <- tau_a*Ia - omega_a*Ra - mu*Ra  # recovered from asymptomatic infection
    
    dBH <- u - b*v_H*BH - (1-b)*v_H*BH # hyper-infectious bacteria
    dBL <- n*BL*(1-(BL/K_san)) + b*v_H*BH - v_L*BL # natural-state bacteria
    
    
    dCDc <- (1-tau_t)*delta_1*theta_12*I1 + (1-tau_t)*delta_2*theta_23*I2 # count compartment: cumulative cholera deaths
    dCInc <- (lam_e+lam_h)*S + (1-l1)^a*(lam_e+lam_h)*V1 + (1-l2)^a*(lam_e+lam_h)*V2 # count compartment: cumulative incidence
    
    # count compartment: cumulative reported deaths
    dCRD <- ((rr_cd * (1-tau_t)*delta_1*theta_12*I1) +    # reported I1 cases --> cholera deaths
               (rr_cd * (1-tau_t)*delta_2*theta_23*I2))  # reported I2 cases --> cholera deaths
    
    # count compartment: cumulative reported incidence
    dCRInc <- ((rr_cd * (1-tau_t)*delta_1*theta_12*I1) +       # reported I1 cases --> cholera deaths 
                 (rr_t * tau_t*theta_12*I1) +                  # reported I1 cases --> treated 
                 (rr_cd * (1-tau_t)*delta_2*theta_23*I2) +     # reported I2 cases --> cholera deaths 
                 (rr_t * tau_t*theta_23*I2))                   # reported I2 cases --> treated
    
    # count comapartment: cumulative severe/moderate incidence
    dCSM <- j*gamma*(E + q*Ev) +  f*gamma*(E + q*Ev)
    
    
    
    # model outputs ----
    list( 
      c(dS, dV1, dV2, dE, dEv, dI1, dI2, dI3, dI4, dIa, dRc, dRcp, dRa, dBH, dBL, dCDc, dCInc, dCRD, dCRInc, dCSM), 
      temperature = temp,
      rainfall = rain,
      n = n, 
      beta_e = beta_e,
      beta_e0 = beta_e0,
      lam_e = lam_e,
      lam_h = lam_h,
      u = u,
      P = P,
      rho1 = rho1,
      rho2 = rho2,
      HYG = HYG,
      SAN = SAN, 
      CWA = CWA,
      hyg_cov = hyg_cov,
      san_cov = san_cov,
      cwa_cov = cwa_cov,
      e_inc = lam_e*S,
      h_inc = lam_h*S,
      K_san = K_san,
      KL = KL)
  })
}


### inital conditions (match to compartment name) ----
initP <- 17700000 # initial population 
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
initS <- initP - initV1 - initV2 - initE - initEv - initI1 - initI2 - initI3 - initI4 - initIa - initRc - initRcp - initRa # calculate intial suscpetible from other conditions

initBH <- 1.48e11
initBL <- 1.92e10

initCDc <- 0
initCInc <- 0
initCRD <- 0
initCRInc <- 0
initCSM <- 0

istate <- c(S = initS, # associate initial consition with the compartment it represents
            V1 = initV1, V2 = initV2,
            E = initE, Ev = initEv,
            I1 = initI1, I2 = initI2, I3 = initI3, I4 = initI4, Ia = initIa,
            Rc = initRc, Rcp = initRcp, Ra = initRa, 
            BH = initBH, BL = initBL,
            CDc = initCDc, CInc = initCInc,
            CRD = initCRD, CRInc = initCRInc, CSM = initCSM)


### parameters ----
parameters <- list(
  mu = 1/(67.8*365), # human natural birth/death rate 
  
  beta_e0 = 2.3e-4, # average environmental infectious contact rate
  beta_h = 4.3e-6, # human infectious contact rate
  
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
  
  tau_t = 0.80, # treatment rate of severe/moderate cases
  
  delta_1 = 0.5, # cholera death rate of untreated severe cases
  delta_2 = 0.1, # cholera death rate of untreated moderate cases
  
  zeta = 1e9, # severe case human bacterial shedding rate
  sigma_2 = 0.1, # relative shedding rate for moderate case
  sigma_3 = 0.01, # relative shedding rate for mild case 
  sigma_4 = 0.0001, # relative shedding rate for asymptomatic case
  
  rr_cd = 0.29, # reporting rate of deaths due to cholera
  rr_t = 1, # reporting rate of cases who were treated
  
  v_H = 5, # hyper-infectious bacteria natural death rate/state-transition rate
  b = 0.1, # proportion of hyper-infectious bacteria secreted into the natural environment
  
  n_0 = 0.33, # mean natural-state bacterial intrinsic proliferation rate
  v_L = 0.33, # natural-state bacteria natural death rate
  KL = 1e7, # concentration of natural-state V. cholerae in water that yields 50% chance of infection
  K = 1e10, # environmental carrying capacity concentration of natural-state bacteria
  
  temp_mean = 21.99, # mean temperature (celsius)
  temp_max = 23.98, # maximum temperature (celsius)
  rain_max = 1067245, # maximum rainfall (mm)
  A_temp = 1.18, # amplitude temperature equation
  A_rain = 319879.08, # amplitude rainfall equation
  phi_temp = 7.16, # phase shift temperature equation
  phi_rain = 8.87, # phase shift rainfall equation
  rain_mean = 414168.02, # mean rainfall (mm)
  
  alpha_n = 1.25, # natural state bacterial intrinsic proliferation sensitivity to temperature
  alpha_b = 2, # environmental infectious contact rate sensitivity to rainfall
  lag_tr = 5, #temperature/rainfall lag
  
  omega_v1 = 1/(2*365), # one-dose loss of vaccine-induced immunity
  omega_v2 = 1/(6.5*365), # 2-dose loss of vaccine-induced immunity
  a = 2, # vaccine dose assortativity term
  l1 = 0.587, # one-dose vaccine efficacy 
  l2 = 0.65, # two-dose vaccine efficacy
  q = 0.45, # vaccine risk reduction of severe+moderate cases
  r = 0.4, # vaccine risk reduction of mild cases
  
  
  # OCV campaign parameters; one index/campaign across all vectors (for vaccination()) (centered at Jan 2023)
  # campaign start
  vax_starts = c(1470-1306, 1470-1289, 1470-759, 1470-610, 1470-590, 1470-432, 1470-397, 1470-277, 1470-247, 1470-17, 1470+105, 1470+190, 1470+226, 1470+281, 1470+299, 1900),
  # two-dose format?
  vax2D = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  # time to administer one-dose
  vax_dur1 = c(3,4,8,6,7,4,9,7,8,7,9,10,7,10,6,10),
  # people vaccinated with one dose
  pvax1 = c(72500, 412767, 1544464, 1167939, 840774, 101245, 1036869, 160340, 265188, 100713, 1909405, 2229941, 1858472, 1522407, 862326, 4406752),
  # time to administer second dose
  vax_dur2 = c(3, 10, 6, 7, NA, 4, 6, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  # interval between first and second dose
  vax_int = c(14, 168, 63, 84, NA, 14, 91, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  # people vaccinated with second dose
  pvax2 = c(72500, 286203, 1515214, 1164639, NA, 101245, 1036859, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  
  
  # WASH coverage improvement campaign parameters; one index/campaign across vectors 
  # for hygiene_coverage()
  hyg_0 = 0.05, # initial coverage
  hyg_starts = c(), # campaign start times
  hyg_ends = c(1001), #3270 # campaign end times
  hyg_goals = c(0.8), # goal coverages 
  
  # for sanitation_coverage()
  san_0 = 0.09, # initial coverage
  san_starts = c(), # campaign start times
  san_ends = c(1001), #3270 # campaign end times
  san_goals = c(0.8), # goal coverages 
  
  # for cleanwater_coverage()
  cwa_0 = 0.5, # intial coverage
  cwa_starts = c(), # campaign start times
  cwa_ends = c(1001), #3270 # campaign end times
  cwa_goals = c(0.9) # goal coverages
)


### timestep ----
time_start <- 1
time_stop <- 365 * 10 # 10 year simulation   
deltat <- 1 # 1 day

tps <- seq(time_start, time_stop, by = deltat) 


### base model output (daily) ----
out <- deSolve::ode(y = istate, times = tps, func = cholera_dpm_vax, parms = parameters) # calculate odes
out_c <- as.data.frame(out)
out_c$I <- out_c$I1+out_c$I2+out_c$I3+out_c$I4+out_c$Ia
out_c$R <- out_c$Rc+out_c$Rcp+out_c$Ra


### relevant model output (monthly) ----
source("scripts/sim_cholera_inc.R") # source function that convert outputs from daily timestep to a data comparable format

dat <- read_excel("data/raw_data/Cholera_Data_Monthly_Mule_1.xlsx") # import cholera case data
CRInc <- out_c$CRInc
CRD <- out_c$CRD

# create data frame containing data and model outputs in comparable formats 
inc_sim.base <- sim_cholera_inc(CRInc = CRInc, CRD = CRD, 
                                c_data = dat$incidence[6:length(dat$incidence)],
                                d_data = dat$death[6:length(dat$death)],
                                years_init = 4, time_stop = time_stop) 

