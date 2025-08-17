### sim_cholera_inc.R
#
# This script hold a function that takes in reported incidence and death data from base model outputs,
# observed incidence and death data from Ethiopia, number of initializations years and a time stop
#
# This script produce a data frame with columns for months, monthly cases incidence, monthly deaths and monthly cfr 
# for the model and the observed data for comparison 
#
###

### define function
sim_cholera_inc <- function(CRInc, CRD, c_data, d_data, years_init, time_stop) {
  months_sim <- seq(1, floor(length(CRInc)/30)) # for a model in days, how many months of data were simulated
  
  # cumulative values for first month on day 30  
  m1_cinc <- CRInc[30] 
  m1_dinc <- CRD[30]
  
  # values for first month first index of results vectors
  calc_cinc <- c(m1_cinc) 
  calc_dinc <- c(m1_dinc)
  
  # for each month simulated calculate the incidence by subtracting the cumulative incidence at the last day of the month - the first day of the month
  for(i in 1:length(months_sim)) { 
    valc <-  CRInc[30*i] - CRInc[30*(i-1)] # incidence
    calc_cinc <- c(calc_cinc,valc)
    
    vald <-  CRD[30*i] - CRD[30*(i-1)] # deaths
    calc_dinc <- c(calc_dinc,vald)
  }
  
  inc_sim <- data.frame(month = months_sim, rcaseinc_sim = calc_cinc, rdeathinc_sim = calc_dinc) # save simulated values to data frame
  
  # add NAs before observed data to match it with a specific starting point in simulated outputs
  # add NAs after the last data point until simulated and obsrved vectors are the same length
  start <- years_init*12+1
  cinc_data <- rep(NA, floor(time_stop/30))
  cinc_data[start:(start+length(c_data)-1)] <- c_data # for incidence
  
  start <- years_init*12+1
  dinc_data <- rep(NA, floor(time_stop/30))
  dinc_data[start:(start+length(d_data)-1)] <- d_data # for deaths

  # save observed values to data frame
  inc_sim$rcaseinc_data <- cinc_data
  inc_sim$rdeathinc_data <- dinc_data
  
  # calculate simulated and observed cfrs and add to data frame
  inc_sim$cfr_sim <- calc_dinc/calc_cinc
  inc_sim$cfr_data <- dinc_data/cinc_data
  
  return(inc_sim)
}
