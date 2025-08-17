### run_scenarioanalysis.R ----
#
# Given the transmission model and base parameters determined from model calibration,
# this code updates relevant parameters depending on a given scenario, runs the transmission model
# with these new parameters and calculates prospective analysis outputs of interest including,
# 1) total case incidence between Jan2025-Jan2028, 2) deaths between Jan2025-Jan2028, and 3) severe 
# +moderate case incidence between Jan2025-2028. These outputs are then saved in a .rds file. 
#
###

### load relevant packages
library(deSolve)


### create data frame to store scenario outputs
scenario_results <- data.frame(scenario = c(), incidence25_28 = c(), deaths25_28 = c(), sevmod_incidence25_28 = c())


### source the base model 
source("scripts/cholera_dpm_vax.R")  # source transmission model for runs, and base parameters to change, given relevant scenario


### define and run scenarios 

## scenario 1 :: the plan (++) :: Full Doses; 2D Preventative; Improve All WASH ----
# update base parameters to fulfill scenario 1
parameters_s1 <- parameters

parameters_s1$vax_starts <- c(1470-1306, 1470-1289, 1470-759, 810, 1170, 1530, 1890, 2250) #
parameters_s1$vax2D <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) #
parameters_s1$vax_dur1 <- c(3, 4, 8, 6, 6, 6, 6, 6) #
parameters_s1$pvax1 <- c(72500, 412767, 1544464, 3031267, 2998115, 3001158, 2972525, 3488851) #
parameters_s1$vax_dur2 <- c(3, 10, 6, 6, 6, 6, 6, 6) #
parameters_s1$vax_int <- c(14, 168, 63, 72, 72, 72, 72, 72) #
parameters_s1$pvax2 <- c(72500, 286203, 1515214, 3031267, 2998115, 3001158, 2972525, 3488851) #

parameters_s1$hyg_starts <- c(1110) #
parameters_s1$hyg_ends <- c(3270) #
parameters_s1$hyg_goals <- c(0.8) #

parameters_s1$san_starts <- c(1110) #
parameters_s1$san_ends <- c(3270) #
parameters_s1$san_goals <- c(0.8) #

parameters_s1$cwa_starts <- c(1110) #
parameters_s1$cwa_ends <- c(3270) #
parameters_s1$cwa_goals <- c(0.9) #

# run the model using scenario 1 parameters
out_s1 <- deSolve::ode(y = istate, times = tps, func = cholera_dpm_vax, parms = parameters_s1) 
out_s1 <- as.data.frame(out_s1)

# calculate outputs of interest 
s1_incidence <- out_s1$CInc[109*30] - out_s1$CInc[73*30] # case incidence between Jan2025 and Jan 2028
s1_deaths <- out_s1$CDc[109*30] - out_s1$CDc[73*30] # deaths between Jan2025 and Jan2028
s1_sevmod <- out_s1$CSM[109*30] - out_s1$CSM[73*30] # severe and moderate case incidence between Jan2025 and Jan2028

# update scenario analysis results data frame with scenario 1 specific outputs
s1_results <- list(scenario = 1, incidence25_28 = s1_incidence, deaths25_28 = s1_deaths, sevmod_incidence25_28 = s1_sevmod)
scenario_results <- rbind(scenario_results, s1_results)


## scenario 2 :: the reality (--) :: Limited Doeses; 2D Preventative; No WASH  ----
# update base parameters to fulfill scenario 2
parameters_s2 <- parameters

parameters_s2$vax_starts <-  c(1470-1306, 1470-1289, 1470-759, 1470-610, 1470-590, 1470-432, 1470-397, 1470-277, 1470-247, 1470-17, 1470+105, 1470+190, 1470+226, 1470+281, 1470+299, 1900, 2250, 2610, 2970) #
parameters_s2$vax2D <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE) #
parameters_s2$vax_dur1 <- c(3,4,8,6,7,4,9,7,8,7,9,10,7,10,6,10, 6, 6, 6) #
parameters_s2$pvax1 <- c(72500, 412767, 1544464, 1167939, 840774, 101245, 1036869, 160340, 265188, 100713, 1909405, 2229941, 1858472, 1522407, 862326, 4406752, 2203376, 2203376, 2203376) #
parameters_s2$vax_dur2 <- c(3, 10, 6, 7, NA, 4, 6, NA, NA, NA, NA, NA, NA, NA, NA, NA, 6, 6, 6) #
parameters_s2$vax_int <- c(14, 168, 63, 84, NA, 14, 91, NA, NA, NA, NA, NA, NA, NA, NA, NA, 72, 72, 72) #
parameters_s2$pvax2 <- c(72500, 286203, 1515214, 1164639, NA, 101245, 1036859, NA, NA, NA, NA, NA, NA, NA, NA, NA, 2203376, 2203376, 2203376) #

parameters_s2$hyg_starts <- numeric(0) #
parameters_s2$hyg_ends <- numeric(0) #
parameters_s2$hyg_goals <- numeric(0) #

parameters_s2$san_starts <- numeric(0) #
parameters_s2$san_ends <- numeric(0) #
parameters_s2$san_goals <- numeric(0) #

parameters_s2$cwa_starts <- numeric(0) #
parameters_s2$cwa_ends <- numeric(0) #
parameters_s2$cwa_goals <- numeric(0) #

# run the model using scenario 2 parameters
out_s2 <- deSolve::ode(y = istate, times = tps, func = cholera_dpm_vax, parms = parameters_s2) 
out_s2 <- as.data.frame(out_s2)

# calculate outputs of interest 
s2_incidence <- out_s2$CInc[109*30] - out_s2$CInc[73*30] # case incidence between Jan2025 and Jan 2028
s2_deaths <- out_s2$CDc[109*30] - out_s2$CDc[73*30]  # deaths between Jan2025 and Jan2028
s2_sevmod <- out_s2$CSM[109*30] - out_s2$CSM[73*30]  # severe and moderate case incidence between Jan2025 and Jan2028

# update scenario analysis results data frame with scenario 2 specific outputs
s2_results <- list(scenario = 2, incidence25_28 = s2_incidence, deaths25_28 = s2_deaths, sevmod_incidence25_28 = s2_sevmod)
scenario_results <- rbind(scenario_results, s2_results)


## scenario 3 :: the reality (--) :: Limited Doses; 1D Reactive; No WASH [BASELINE]  ----
# update base parameters to fulfill scenario 3
parameters_s3<- parameters

parameters_s3$vax_starts <-  c(1470-1306, 1470-1289, 1470-759, 1470-610, 1470-590, 1470-432, 1470-397, 1470-277, 1470-247, 1470-17, 1470+105, 1470+190, 1470+226, 1470+281, 1470+299, 1900, 2460, 2820, 3180) #
parameters_s3$vax2D <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) #
parameters_s3$vax_dur1 <- c(3,4,8,6,7,4,9,7,8,7,9,10,7,10,6,10, 6, 6, 6) #
parameters_s3$pvax1 <- c(72500, 412767, 1544464, 1167939, 840774, 101245, 1036869, 160340, 265188, 100713, 1909405, 2229941, 1858472, 1522407, 862326, 4406752, 4406752, 4406752, 4406752) #
parameters_s3$vax_dur2 <- c(3, 10, 6, 7, NA, 4, 6, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA) #
parameters_s3$vax_int <- c(14, 168, 63, 84, NA, 14, 91, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA) #
parameters_s3$pvax2 <- c(72500, 286203, 1515214, 1164639, NA, 101245, 1036859, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA) #

parameters_s3$hyg_starts <- numeric(0)
parameters_s3$hyg_ends <- numeric(0)
parameters_s3$hyg_goals <- numeric(0)

parameters_s3$san_starts <- numeric(0)
parameters_s3$san_ends <- numeric(0)
parameters_s3$san_goals <- numeric(0)

parameters_s3$cwa_starts <- numeric(0)
parameters_s3$cwa_ends <- numeric(0)
parameters_s3$cwa_goals <- numeric(0)

# run the model using scenario 3 parameters
out_s3 <- deSolve::ode(y = istate, times = tps, func = cholera_dpm_vax, parms = parameters_s3) 
out_s3 <- as.data.frame(out_s3)

# calculate outputs of interest 
s3_incidence <- out_s3$CInc[109*30] - out_s3$CInc[73*30] # case incidence between Jan2025 and Jan 2028
s3_deaths <- out_s3$CDc[109*30] - out_s3$CDc[73*30] # deaths between Jan2025 and Jan2028
s3_sevmod <- out_s3$CSM[109*30] - out_s3$CSM[73*30] # severe and moderate case incidence between Jan2025 and Jan2028

# update scenario analysis results data frame with scenario 3 specific outputs
s3_results <- list(scenario = 3, incidence25_28 = s3_incidence, deaths25_28 = s3_deaths, sevmod_incidence25_28 = s3_sevmod)
scenario_results <- rbind(scenario_results, s3_results)



## scenario 4 :: the possibility (-+) :: Limited Doses; 1D Reactive; Improve HYG  ----
# update base parameters to fulfill scenario 4
parameters_s4<- parameters

parameters_s4$vax_starts <-  c(1470-1306, 1470-1289, 1470-759, 1470-610, 1470-590, 1470-432, 1470-397, 1470-277, 1470-247, 1470-17, 1470+105, 1470+190, 1470+226, 1470+281, 1470+299, 1900, 2460, 2820, 3180)
parameters_s4$vax2D <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
parameters_s4$vax_dur1 <- c(3,4,8,6,7,4,9,7,8,7,9,10,7,10,6,10, 6, 6, 6)
parameters_s4$pvax1 <- c(72500, 412767, 1544464, 1167939, 840774, 101245, 1036869, 160340, 265188, 100713, 1909405, 2229941, 1858472, 1522407, 862326, 4406752, 4406752, 4406752, 4406752) #
parameters_s4$vax_dur2 <- c(3, 10, 6, 7, NA, 4, 6, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s4$vax_int <- c(14, 168, 63, 84, NA, 14, 91, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s4$pvax2 <- c(72500, 286203, 1515214, 1164639, NA, 101245, 1036859, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

parameters_s4$hyg_starts <- c(2190) #
parameters_s4$hyg_ends <- c(3270) #
parameters_s4$hyg_goals <- c(0.8) #

parameters_s4$san_starts <- numeric(0)
parameters_s4$san_ends <- numeric(0)
parameters_s4$san_goals <- numeric(0)

parameters_s4$cwa_starts <- numeric(0)
parameters_s4$cwa_ends <- numeric(0)
parameters_s4$cwa_goals <- numeric(0)

# run the model using scenario 4 parameters
out_s4 <- deSolve::ode(y = istate, times = tps, func = cholera_dpm_vax, parms = parameters_s4) 
out_s4 <- as.data.frame(out_s4)

# calculate outputs of interest 
s4_incidence <- out_s4$CInc[109*30] - out_s4$CInc[73*30] # case incidence between Jan2025 and Jan 2028
s4_deaths <- out_s4$CDc[109*30] - out_s4$CDc[73*30] # deaths between Jan2025 and Jan2028
s4_sevmod <- out_s4$CSM[109*30] - out_s4$CSM[73*30] # severe and moderate case incidence between Jan2025 and Jan2028

# update scenario analysis results data frame with scenario 4 specific outputs
s4_results <- list(scenario = 4, incidence25_28 = s4_incidence, deaths25_28 = s4_deaths, sevmod_incidence25_28 = s4_sevmod)
scenario_results <- rbind(scenario_results, s4_results)



## scenario 5 :: the possibility (-+) :: Limited Doses; 1D Reactive; Improve SAN  ----
# update base parameters to fulfill scenario 5
parameters_s5<- parameters

parameters_s5$vax_starts <-  c(1470-1306, 1470-1289, 1470-759, 1470-610, 1470-590, 1470-432, 1470-397, 1470-277, 1470-247, 1470-17, 1470+105, 1470+190, 1470+226, 1470+281, 1470+299, 1900, 2460, 2820, 3180)
parameters_s5$vax2D <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
parameters_s5$vax_dur1 <- c(3,4,8,6,7,4,9,7,8,7,9,10,7,10,6,10, 6, 6, 6)
parameters_s5$pvax1 <- c(72500, 412767, 1544464, 1167939, 840774, 101245, 1036869, 160340, 265188, 100713, 1909405, 2229941, 1858472, 1522407, 862326, 4406752, 4406752, 4406752, 4406752)
parameters_s5$vax_dur2 <- c(3, 10, 6, 7, NA, 4, 6, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s5$vax_int <- c(14, 168, 63, 84, NA, 14, 91, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s5$pvax2 <- c(72500, 286203, 1515214, 1164639, NA, 101245, 1036859, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

parameters_s5$hyg_starts <- numeric(0)
parameters_s5$hyg_ends <- numeric(0)
parameters_s5$hyg_goals <- numeric(0)

parameters_s5$san_starts <- c(2190) #
parameters_s5$san_ends <- c(3270) #
parameters_s5$san_goals <- c(0.8) #

parameters_s5$cwa_starts <- numeric(0)
parameters_s5$cwa_ends <- numeric(0)
parameters_s5$cwa_goals <- numeric(0)

# run the model using scenario 5 parameters
out_s5 <- deSolve::ode(y = istate, times = tps, func = cholera_dpm_vax, parms = parameters_s5) 
out_s5 <- as.data.frame(out_s5)

# calculate outputs of interest 
s5_incidence <- out_s5$CInc[109*30] - out_s5$CInc[73*30] # case incidence between Jan2025 and Jan 2028
s5_deaths <- out_s5$CDc[109*30] - out_s5$CDc[73*30] # deaths between Jan2025 and Jan2028
s5_sevmod <- out_s5$CSM[109*30] - out_s5$CSM[73*30] # severe and moderate case incidence between Jan2025 and Jan2028

# update scenario analysis results data frame with scenario 5 specific outputs
s5_results <- list(scenario = 5, incidence25_28 = s5_incidence, deaths25_28 = s5_deaths, sevmod_incidence25_28 = s5_sevmod)
scenario_results <- rbind(scenario_results, s5_results)



## scenario 6 :: the possibility (-+) :: Limited Doses; 1D Reactive; Improve CWA ----
# update base parameters to fulfill scenario 6
parameters_s6<- parameters

parameters_s6$vax_starts <-  c(1470-1306, 1470-1289, 1470-759, 1470-610, 1470-590, 1470-432, 1470-397, 1470-277, 1470-247, 1470-17, 1470+105, 1470+190, 1470+226, 1470+281, 1470+299, 1900, 2460, 2820, 3180)
parameters_s6$vax2D <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
parameters_s6$vax_dur1 <- c(3,4,8,6,7,4,9,7,8,7,9,10,7,10,6,10, 6, 6, 6)
parameters_s6$pvax1 <- c(72500, 412767, 1544464, 1167939, 840774, 101245, 1036869, 160340, 265188, 100713, 1909405, 2229941, 1858472, 1522407, 862326, 4406752, 4406752, 4406752, 4406752)
parameters_s6$vax_dur2 <- c(3, 10, 6, 7, NA, 4, 6, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s6$vax_int <- c(14, 168, 63, 84, NA, 14, 91, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s6$pvax2 <- c(72500, 286203, 1515214, 1164639, NA, 101245, 1036859, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

parameters_s6$hyg_starts <- numeric(0)
parameters_s6$hyg_ends <- numeric(0)
parameters_s6$hyg_goals <- numeric(0)

parameters_s6$san_starts <- numeric(0)
parameters_s6$san_ends <- numeric(0)
parameters_s6$san_goals <- numeric(0)

parameters_s6$cwa_starts <- c(2190) #
parameters_s6$cwa_ends <- c(3270) #
parameters_s6$cwa_goals <- c(0.9) #

# run the model using scenario 6 parameters
out_s6 <- deSolve::ode(y = istate, times = tps, func = cholera_dpm_vax, parms = parameters_s6) 
out_s6 <- as.data.frame(out_s6)

# calculate outputs of interest 
s6_incidence <- out_s6$CInc[109*30] - out_s6$CInc[73*30] # case incidence between Jan2025 and Jan 2028
s6_deaths <- out_s6$CDc[109*30] - out_s6$CDc[73*30] # deaths between Jan2025 and Jan2028
s6_sevmod <- out_s6$CSM[109*30] - out_s6$CSM[73*30] # severe and moderate case incidence between Jan2025 and Jan2028

# update scenario analysis results data frame with scenario 6 specific outputs
s6_results <- list(scenario = 6, incidence25_28 = s6_incidence, deaths25_28 = s6_deaths, sevmod_incidence25_28 = s6_sevmod)
scenario_results <- rbind(scenario_results, s6_results)



## scenario 7 :: the possibility (-+) :: Limited Doses; 1D Reactive; Improve HYG/SAN ----
# update base parameters to fulfill scenario 7
parameters_s7<- parameters

parameters_s7$vax_starts <-  c(1470-1306, 1470-1289, 1470-759, 1470-610, 1470-590, 1470-432, 1470-397, 1470-277, 1470-247, 1470-17, 1470+105, 1470+190, 1470+226, 1470+281, 1470+299, 1900, 2460, 2820, 3180)
parameters_s7$vax2D <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
parameters_s7$vax_dur1 <- c(3,4,8,6,7,4,9,7,8,7,9,10,7,10,6,10, 6, 6, 6)
parameters_s7$pvax1 <- c(72500, 412767, 1544464, 1167939, 840774, 101245, 1036869, 160340, 265188, 100713, 1909405, 2229941, 1858472, 1522407, 862326, 4406752, 4406752, 4406752, 4406752)
parameters_s7$vax_dur2 <- c(3, 10, 6, 7, NA, 4, 6, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s7$vax_int <- c(14, 168, 63, 84, NA, 14, 91, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s7$pvax2 <- c(72500, 286203, 1515214, 1164639, NA, 101245, 1036859, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

parameters_s7$hyg_starts <- c(2190) #
parameters_s7$hyg_ends <- c(3270) #
parameters_s7$hyg_goals <- c(0.8) #

parameters_s7$san_starts <- c(2190) #
parameters_s7$san_ends <- c(3270) #
parameters_s7$san_goals <- c(0.8) #

parameters_s7$cwa_starts <- numeric(0)
parameters_s7$cwa_ends <- numeric(0)
parameters_s7$cwa_goals <- numeric(0)

# run the model using scenario 7 parameters
out_s7 <- deSolve::ode(y = istate, times = tps, func = cholera_dpm_vax, parms = parameters_s7) 
out_s7 <- as.data.frame(out_s7)

# calculate outputs of interest 
s7_incidence <- out_s7$CInc[109*30] - out_s7$CInc[73*30] # case incidence between Jan2025 and Jan 2028
s7_deaths <- out_s7$CDc[109*30] - out_s7$CDc[73*30] # deaths between Jan2025 and Jan2028
s7_sevmod <- out_s7$CSM[109*30] - out_s7$CSM[73*30] # severe and moderate case incidence between Jan2025 and Jan2028

# update scenario analysis results data frame with scenario 7 specific outputs
s7_results <- list(scenario = 7, incidence25_28 = s7_incidence, deaths25_28 = s7_deaths, sevmod_incidence25_28 = s7_sevmod)
scenario_results <- rbind(scenario_results, s7_results)



## scenario 8 :: the possibility (-+) :: Limited Doses; 1D Reactive; Improve HYG/CWA ----
# update base parameters to fulfill scenario 8
parameters_s8 <- parameters

parameters_s8$vax_starts <-  c(1470-1306, 1470-1289, 1470-759, 1470-610, 1470-590, 1470-432, 1470-397, 1470-277, 1470-247, 1470-17, 1470+105, 1470+190, 1470+226, 1470+281, 1470+299, 1900, 2460, 2820, 3180)
parameters_s8$vax2D <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
parameters_s8$vax_dur1 <- c(3,4,8,6,7,4,9,7,8,7,9,10,7,10,6,10, 6, 6, 6)
parameters_s8$pvax1 <- c(72500, 412767, 1544464, 1167939, 840774, 101245, 1036869, 160340, 265188, 100713, 1909405, 2229941, 1858472, 1522407, 862326, 4406752, 4406752, 4406752, 4406752)
parameters_s8$vax_dur2 <- c(3, 10, 6, 7, NA, 4, 6, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s8$vax_int <- c(14, 168, 63, 84, NA, 14, 91, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s8$pvax2 <- c(72500, 286203, 1515214, 1164639, NA, 101245, 1036859, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

parameters_s8$hyg_starts <- c(2190)
parameters_s8$hyg_ends <- c(3270)
parameters_s8$hyg_goals <- c(0.8)

parameters_s8$san_starts <- numeric(0) #
parameters_s8$san_ends <- numeric(0) #
parameters_s8$san_goals <- numeric(0) #

parameters_s8$cwa_starts <- c(2190) #
parameters_s8$cwa_ends <- c(3270) #
parameters_s8$cwa_goals <- c(0.9) #

# run the model using scenario 8 parameters
out_s8 <- deSolve::ode(y = istate, times = tps, func = cholera_dpm_vax, parms = parameters_s8) 
out_s8 <- as.data.frame(out_s8)

# calculate outputs of interest 
s8_incidence <- out_s8$CInc[109*30] - out_s8$CInc[73*30] # case incidence between Jan2025 and Jan 2028
s8_deaths <- out_s8$CDc[109*30] - out_s8$CDc[73*30] # deaths between Jan2025 and Jan2028
s8_sevmod <- out_s8$CSM[109*30] - out_s8$CSM[73*30] # severe and moderate case incidence between Jan2025 and Jan2028

# update scenario analysis results data frame with scenario 8 specific outputs
s8_results <- list(scenario = 8, incidence25_28 = s8_incidence, deaths25_28 = s8_deaths, sevmod_incidence25_28 = s8_sevmod)
scenario_results <- rbind(scenario_results, s8_results)



## scenario 9 :: the possibility (-+) :: Limited Doses; 1D Reactive; Improve SAN/CWA ----
# update base parameters to fulfill scenario 9
parameters_s9 <- parameters

parameters_s9$vax_starts <-  c(1470-1306, 1470-1289, 1470-759, 1470-610, 1470-590, 1470-432, 1470-397, 1470-277, 1470-247, 1470-17, 1470+105, 1470+190, 1470+226, 1470+281, 1470+299, 1900, 2460, 2820, 3180)
parameters_s9$vax2D <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
parameters_s9$vax_dur1 <- c(3,4,8,6,7,4,9,7,8,7,9,10,7,10,6,10, 6, 6, 6)
parameters_s9$pvax1 <- c(72500, 412767, 1544464, 1167939, 840774, 101245, 1036869, 160340, 265188, 100713, 1909405, 2229941, 1858472, 1522407, 862326, 4406752, 4406752, 4406752, 4406752)
parameters_s9$vax_dur2 <- c(3, 10, 6, 7, NA, 4, 6, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s9$vax_int <- c(14, 168, 63, 84, NA, 14, 91, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s9$pvax2 <- c(72500, 286203, 1515214, 1164639, NA, 101245, 1036859, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

parameters_s9$hyg_starts <- numeric(0) #
parameters_s9$hyg_ends <- numeric(0) #
parameters_s9$hyg_goals <- numeric(0) #

parameters_s9$san_starts <- c(2190) 
parameters_s9$san_ends <- c(3270) 
parameters_s9$san_goals <- c(0.8) #

parameters_s9$cwa_starts <- c(2190)
parameters_s9$cwa_ends <- c(3270)
parameters_s9$cwa_goals <- c(0.9)

# run the model using scenario 9 parameters
out_s9 <- deSolve::ode(y = istate, times = tps, func = cholera_dpm_vax, parms = parameters_s9) 
out_s9 <- as.data.frame(out_s9)

# calculate outputs of interest 
s9_incidence <- out_s9$CInc[109*30] - out_s9$CInc[73*30] # case incidence between Jan2025 and Jan 2028
s9_deaths <- out_s9$CDc[109*30] - out_s9$CDc[73*30] # deaths between Jan2025 and Jan2028
s9_sevmod <- out_s9$CSM[109*30] - out_s9$CSM[73*30] # severe and moderate case incidence between Jan2025 and Jan2028

# update scenario analysis results data frame with scenario 9 specific outputs
s9_results <- list(scenario = 9, incidence25_28 = s9_incidence, deaths25_28 = s9_deaths, sevmod_incidence25_28 = s9_sevmod)
scenario_results <- rbind(scenario_results, s9_results)



## scenario 10 :: the possibility (-+) :: Limited Doses; 1D Reactive; Improve All WASH ----
# update base parameters to fulfill scenario 10
parameters_s10 <- parameters

parameters_s10$vax_starts <-  c(1470-1306, 1470-1289, 1470-759, 1470-610, 1470-590, 1470-432, 1470-397, 1470-277, 1470-247, 1470-17, 1470+105, 1470+190, 1470+226, 1470+281, 1470+299, 1900, 2460, 2820, 3180)
parameters_s10$vax2D <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
parameters_s10$vax_dur1 <- c(3,4,8,6,7,4,9,7,8,7,9,10,7,10,6,10, 6, 6, 6)
parameters_s10$pvax1 <- c(72500, 412767, 1544464, 1167939, 840774, 101245, 1036869, 160340, 265188, 100713, 1909405, 2229941, 1858472, 1522407, 862326, 4406752, 4406752, 4406752, 4406752)
parameters_s10$vax_dur2 <- c(3, 10, 6, 7, NA, 4, 6, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s10$vax_int <- c(14, 168, 63, 84, NA, 14, 91, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s10$pvax2 <- c(72500, 286203, 1515214, 1164639, NA, 101245, 1036859, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

parameters_s10$hyg_starts <- c(2190) #
parameters_s10$hyg_ends <- c(3270) #
parameters_s10$hyg_goals <- c(0.8) #

parameters_s10$san_starts <- c(2190)
parameters_s10$san_ends <- c(3270)
parameters_s10$san_goals <- c(0.8)

parameters_s10$cwa_starts <- c(2190)
parameters_s10$cwa_ends <- c(3270)
parameters_s10$cwa_goals <- c(0.9)

# run the model using scenario 10 parameters
out_s10 <- deSolve::ode(y = istate, times = tps, func = cholera_dpm_vax, parms = parameters_s10) 
out_s10 <- as.data.frame(out_s10)

# calculate outputs of interest 
s10_incidence <- out_s10$CInc[109*30] - out_s10$CInc[73*30] # case incidence between Jan2025 and Jan 2028
s10_deaths <- out_s10$CDc[109*30] - out_s10$CDc[73*30] # deaths between Jan2025 and Jan2028
s10_sevmod <- out_s10$CSM[109*30] - out_s10$CSM[73*30] # severe and moderate case incidence between Jan2025 and Jan2028

# update scenario analysis results data frame with scenario 10 specific outputs
s10_results <- list(scenario = 10, incidence25_28 = s10_incidence, deaths25_28 = s10_deaths, sevmod_incidence25_28 = s10_sevmod)
scenario_results <- rbind(scenario_results, s10_results)



## scenario 11 :: the possibility (-+/-) :: Limited Doses; 1D Reactive; Partial Improvement of all Wash----
# update base parameters to fulfill scenario 10
parameters_s11 <- parameters

parameters_s11$vax_starts <-  c(1470-1306, 1470-1289, 1470-759, 1470-610, 1470-590, 1470-432, 1470-397, 1470-277, 1470-247, 1470-17, 1470+105, 1470+190, 1470+226, 1470+281, 1470+299, 1900, 2460, 2820, 3180)
parameters_s11$vax2D <- c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
parameters_s11$vax_dur1 <- c(3,4,8,6,7,4,9,7,8,7,9,10,7,10,6,10, 6, 6, 6)
parameters_s11$pvax1 <- c(72500, 412767, 1544464, 1167939, 840774, 101245, 1036869, 160340, 265188, 100713, 1909405, 2229941, 1858472, 1522407, 862326, 4406752, 4406752, 4406752, 4406752)
parameters_s11$vax_dur2 <- c(3, 10, 6, 7, NA, 4, 6, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s11$vax_int <- c(14, 168, 63, 84, NA, 14, 91, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
parameters_s11$pvax2 <- c(72500, 286203, 1515214, 1164639, NA, 101245, 1036859, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

parameters_s11$hyg_starts <- c(2190) 
parameters_s11$hyg_ends <- c(3270)
parameters_s11$hyg_goals <- c(0.425) #

parameters_s11$san_starts <- c(2190)
parameters_s11$san_ends <- c(3270)
parameters_s11$san_goals <- c(0.445) #

parameters_s11$cwa_starts <- c(2190)
parameters_s11$cwa_ends <- c(3270)
parameters_s11$cwa_goals <- c(0.7) #

# run the model using scenario 10 parameters
out_s11 <- deSolve::ode(y = istate, times = tps, func = cholera_dpm_vax, parms = parameters_s11) 
out_s11 <- as.data.frame(out_s11)

# calculate outputs of interest 
s11_incidence <- out_s11$CInc[109*30] - out_s11$CInc[73*30] # case incidence between Jan2025 and Jan 2028
s11_deaths <- out_s11$CDc[109*30] - out_s11$CDc[73*30] # deaths between Jan2025 and Jan2028
s11_sevmod <- out_s11$CSM[109*30] - out_s11$CSM[73*30] # severe and moderate case incidence between Jan2025 and Jan2028

# update scenario analysis results data frame with scenario 11 specific outputs
s11_results <- list(scenario = 11, incidence25_28 = s11_incidence, deaths25_28 = s11_deaths, sevmod_incidence25_28 = s11_sevmod)
scenario_results <- rbind(scenario_results, s11_results)


### save scenario outputs of interested as an .RDS file 
saveRDS(scenario_results, file = "outputs/scenarioanalysis_results.rds")

### save all scenario runs base output as an .RDS file
saveRDS(out_s1, file = "outputs/scenarioanalysis_outs1.rds")
saveRDS(out_s2, file = "outputs/scenarioanalysis_outs2.rds")
saveRDS(out_s3, file = "outputs/scenarioanalysis_outs3.rds")
saveRDS(out_s4, file = "outputs/scenarioanalysis_outs4.rds")
saveRDS(out_s5, file = "outputs/scenarioanalysis_outs5.rds")
saveRDS(out_s6, file = "outputs/scenarioanalysis_outs6.rds")
saveRDS(out_s7, file = "outputs/scenarioanalysis_outs7.rds")
saveRDS(out_s8, file = "outputs/scenarioanalysis_outs8.rds")
saveRDS(out_s9, file = "outputs/scenarioanalysis_outs9.rds")
saveRDS(out_s10, file = "outputs/scenarioanalysis_outs10.rds")
saveRDS(out_s11, file = "outputs/scenarioanalysis_outs11.rds")

