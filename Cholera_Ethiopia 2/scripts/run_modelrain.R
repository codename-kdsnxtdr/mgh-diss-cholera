### run_modelrain.R
#
# This script takes in relevant data and formats it for use by modelrain.stan to
# execute parameter estimation using Bayesian MCMC methods. Parameters being estimated
# were used in a sine function representative of rainfall over time fitted to national
# Ethiopian rainfall data.
#
# This script outputs the numerical posterior outputs of the Stan model for parameters of
# interest and plots that can be used to validate the sufficiency of these estimates including, 
# posterior distributions across chains, trace plots and model outputs plotted against observed data.
#
###


### load relevant packages ----
library(rstan)
library(bayesplot)
library(ggplot2)
library(esquisse)
library(tibble)
library(tidyr)
library(gridExtra)
library(dplyr)


### run the Stan model ----
## set Stan options
rstan_options(auto_write = TRUE) # prevent recompilation to C++ to save time in subsequent runs
options(mc.cores = parallel::detectCores()) # enable parallel execution of  MCMC chains to reduce computation time

## format raw data for recognition by Stan model
dat <- as.data.frame(read.csv("data/raw_data/Malaria_data_full_new.csv")) # load data frame containing monthly rain data
rain.vec <- dat[44:67,4] # define relevant data: rainfall data only between August 2022 and June 2024
rain.vec <- rain.vec/100000 # convert rainfall in mm to a smaller value to prevent confusion/model complexity and simpler estimation of priors 

## define fixed parameters
n_months <- length(rain.vec) # number of data points
ts <- seq(1, n_months, by = 1) # define timesteps, equal to number of data points

## compile parameters for Stan data block
data_rain <- list(
  n_months = n_months,
  ts = ts,
  rain = rain.vec
)

## define MCMC iterations and number of chains
niter <- 4000
nwarmup <- 2000
nchains <- 4

## run Stan model
modelrain_samples <- stan(
  "scripts/modelrain.stan", # call Stan script
  iter = niter,
  warmup = nwarmup,
  data = data_rain,
  chains = nchains,
  cores = nchains,
  seed = 31415, # set seed for reproducibility 
  refresh = 100
)


### produce outputs from Stan model run ----

time <- seq(1, 24, by=1) # define time for plotting raw data 
plot(x = time, y = rain.vec*100000) # plot of raw data used for fitting, converted back to mm

pars <- c('A_Ra', 'phi_Ra', 'Ra_mean', 'sigma') # define sampled/transformed parameters of interest

## numerical outputs 
print(modelrain_samples, pars = pars) # print model output :: including...
  # estimated posterior distribution of paramters of interest (mean, with 95% credible interval)
  # standard deviations 
  # n_eff and Rhat metrics to validate convergence
  

## visual outputs
# 1. posterior densities
stan_dens(modelrain_samples, pars=pars, separate_chains = TRUE) # overlapping distributions of chains suggests convergence 

# 2. trace plots
traceplot(modelrain_samples, pars=pars) # dense, uniform oscillations around a mean following a burn-in period suggest convergence

# 3. observed data vs. predicted temperature + 95% credible interval 
# create a data frame that has the pred_rain distribution from the generated quantities block for each time point as a row, and
# its key features including the mean and 95% credible values and the data at that time point as columns 
rain_pred <- cbind(as.data.frame(summary(
  modelrain_samples, pars = "pred_rain", probs = c(0.05, 0.5, 0.95))$summary), time, dat$rainfall[44:67])
colnames(rain_pred) <- make.names(colnames(rain_pred))

# using the previous data frame, plot the model predicted rainfall for each time point and their credible intervals against the observed data
ggplot(rain_pred, mapping = aes(x = time)) +
  geom_ribbon(aes(ymin = X5., ymax = X95., fill = "95% Credible Interval"), alpha = 0.35) +
  geom_line(mapping = aes(y = X50., color = "Predicted"), size = 0.75) +
  geom_point(mapping = aes(y = dat$rainfall[44:67], color = "Observed"), size = 2) +
  labs(title = "Observed and Predicted Rainfall vs. Time", 
       x = "Time (month)", 
       y = "Rainfall in Millimeters",
       fill = "Legend",
       color = "Legend") +
  scale_fill_manual(values = c("95% Credible Interval" = "lightblue4")) +
  scale_color_manual(values = c("Predicted" = "darkblue", "Observed" = "black")) +
  theme(plot.title = element_text(hjust = 0.5))



### misc rainfall plot ----
## plotsS21B ----
### plot of rainfall data from Ethiopia (January2019-December2024)
##### highlighting data point exclusion of points not between August 2022 and June 2024
#
#
# format data
dat_ex_rain <- as.data.frame(cbind(month = dat[,1], rain = dat[,4]))

# group the data to prevent excluded points on either side of included data from connecting to each other via a line 
dat_ex_rain <- dat_ex_rain %>% # group the data based on whether in was included or excluded (before or after period of inclusion)
  mutate(
    Group = case_when( # define what time points are included or not 
      month < 44 ~ "rain_Excluded_1",
      month >= 44 & month <= 67 ~ "rain_Included",
      month > 67 ~ "rain_Excluded_2"
    )
  )

# ensure that the legend still only has two entries, even though data has three groups
dat_ex_rain <- dat_ex_rain %>% # create legend groups based on inclusion vs. exclusion, ignoring what type of exclusion it is
  mutate(
    LegendGroup = ifelse(grepl("Excluded", Group), "Rainfall (Excluded)", "Rainfall (Included)")
  )

group_colors_3 <- c( # color the plots similarly based on inclusion or exclusion
  "rain_Excluded_1" = "#BBBBBB",
  "rain_Included" = "blue3",
  "rain_Excluded_2" = "#BBBBBB"
)

group_colors_2 <- c( # color the legend based on inclusion or exclusion 
  "Rainfall (Excluded)" = "#BBBBBB",
  "Rainfall (Included)" = "blue3"
)

# produce the plot showing all data points over time grouped based on inclusion vs exclusion 
rain_plot <- ggplot(dat_ex_rain, aes(x = month, y = rain)) +
  # annotate boundaries of data inclusion
  geom_vline(xintercept = 44, linewidth = 1.5, colour = "#999999", linetype = "dotted") + 
  annotate("text", x = 43.8, y = max(dat_ex_rain$rain) * 0.95, label = "August 2022", 
           angle = 90, vjust = -0.8, size = 4, color = "#999999", fontface = "bold") +
  geom_vline(xintercept = 67, linewidth = 1.5, colour = "#999999", linetype = "dotted") +
  annotate("text", x = 67.8, y = max(dat_ex_rain$rain) * 0.95, label = "June 2024", 
           angle = 90, vjust = -0.8, size = 4, color = "#999999", fontface = "bold") +
  # plot all three groups as points connected by lines
  geom_line(data = filter(dat_ex_rain, Group == "rain_Excluded_1"),
            aes(color = Group), size = 0.5) +
  geom_line(data = filter(dat_ex_rain, Group == "rain_Included"),
            aes(color = Group), size = 0.5) +
  geom_line(data = filter(dat_ex_rain, Group == "rain_Excluded_2"),
            aes(color = Group), size = 0.5) +
  geom_point(aes(color = Group), size = 3L) +
  # reference line colors
  scale_color_manual(
    values = group_colors_3,
    guide = "none"
  ) +
  # reference legend colors
  geom_point(aes(fill = LegendGroup), shape = 21, size = 4, color = NA, show.legend = TRUE) +
  scale_fill_manual(
    values = group_colors_2,
    name = "",
    labels = names(group_colors_2)
  ) +
  # define axis names 
  labs(
    x = "Month",
    y = "Millimeters"
  ) +
  # visually format plot
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = c(0.4, 0.90),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

rain_plot






