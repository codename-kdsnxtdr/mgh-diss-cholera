### run_modeltemp.R
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

### load relevant packages ----
library(rstan)
library(bayesplot)
library(ggplot2)
library(esquisse)
library(tibble)
library(tidyr)
library(gridExtra)


### run the Stan model ----
## set Stan options
rstan_options(auto_write = TRUE) # prevent recompilation to C++ to save time in subsequent runs
options(mc.cores = parallel::detectCores()) # enable parallel execution of  MCMC chains to reduce computation time

## format raw data for recognition by Stan model
dat <- as.data.frame(read.csv("data/raw_data/Malaria_data_full_new.csv")) # load data frame containing monthly temp data
temp.vec <- dat[44:67,3] # define relevant data: temp data only between August 2022 and June 2024

## define fixed parameters
n_months <- length(temp.vec) # number of data points
ts <- seq(1, n_months, by = 1) # define timesteps, equal to number of data points

## compile parameters for Stan data block
data_temp <- list(
  n_months = n_months,
  ts = ts,
  temps = temp.vec
)

## define MCMC iterations and number of chains
niter <- 4000
nwarmup <- 2000
nchains <- 4

## run Stan model
modeltemp_samples <- stan(
  "scripts/modeltemp.stan", # call Stan script
  iter = niter,
  warmup = nwarmup,
  data = data_temp,
  chains = nchains,
  cores = nchains,
  seed = 31415, # set seed for reproducibility 
  refresh = 100
)



### produce outputs from Stan model run ----

time <- seq(1, 24, by=1) # define time for plotting raw data 
plot(x = time, y = temp.vec) # plot of raw data used for fitting 

pars <- c('A_Te', 'phi_Te', 'Te_mean', 'sigma') # define sampled/transformed parameters of interest


## numerical outputs 
print(modeltemp_samples, pars = pars) # print model output :: including...
  # estimated posterior distribution of paramters of interest (mean, with 95% credible interval)
  # standard deviations 
  # n_eff and Rhat metrics to validate convergence


## visual outputs
# 1. posterior densities
stan_dens(modeltemp_samples, par = pars, separate_chains = TRUE) # overlapping distributions of chains suggests convergence 


# 2. trace plots
traceplot(modeltemp_samples, pars=pars) # dense, uniform oscillations around a mean following a burn-in period suggest convergence


# 3. observed data vs. predicted temperature + 95% credible interval 
# create a data frame that has the pred_temp distribution from the generated quantities block for each time point as a row, and
# its key features including the mean and 95% credible values and the data at that time point as columns 
temp_pred <- cbind(as.data.frame(summary(
  modeltemp_samples, pars = "pred_temp", probs = c(0.05, 0.5, 0.95))$summary), time, dat$temprature[44:67])
colnames(temp_pred) <- make.names(colnames(temp_pred))

# using the previous data frame, plot the model predicted temperatures for each time point ad their credible intervals against the observed data
ggplot(temp_pred, mapping = aes(x = time)) +
  geom_ribbon(aes(ymin = X5., ymax = X95., fill = "95% Credible Interval"), alpha = 0.35) +
  geom_line(mapping = aes(y = X50., color = "Predicted"), size = 0.75) +
  geom_point(mapping = aes(y = dat$temprature[44:67], color = "Observed"), size = 2) +
  labs(title = "Observed and Predicted Temperature vs. Time", 
       x = "Time (month)", 
       y = "Temperature in Celcius",
       fill = "Legend",
       color = "Legend") +
  scale_fill_manual(values = c("95% Credible Interval" = "lightblue4")) +
  scale_color_manual(values = c("Predicted" = "darkblue", "Observed" = "black")) +
  theme(plot.title = element_text(hjust = 0.5))




### misc temperature plot ----
## plotS21B ----
### plot of temperature data from Ethiopia (January2019-December2024)
##### highlighting data point exclusion of points not between August 2022 and June 2024
#
#
# format data
dat_ex_temp <- as.data.frame(cbind(month = dat[,1], temp = dat[,3]))

# group the data to prevent excluded points on either side of included data from connecting to each other via a line 
dat_ex_temp <- dat_ex_temp %>%  # group the data based on whether in was included or excluded (before or after period of inclusion)
  mutate(
    Group = case_when(# define what time points are included or not 
      month < 44 ~ "temp_Excluded_1",
      month >= 44 & month <= 67 ~ "temp_Included",
      month > 67 ~ "temp_Excluded_2"
    )
  )

# ensure that the legend still only has two entries, even though data has three groups
dat_ex_temp <- dat_ex_temp %>% # create legend groups based on inclusion vs. exclusion, ignoring what type of exclusion it is
  mutate(
    LegendGroup = ifelse(grepl("Excluded", Group), "Temperature (Excluded)", "Temperature (Included)")
  )

group_colors_3 <- c( # color the plots similarly based on inclusion or exclusion
  "temp_Excluded_1" = "#BBBBBB",
  "temp_Included" = "orange2",
  "temp_Excluded_2" = "#BBBBBB"
)

group_colors_2 <- c(  # color the legend based on inclusion or exclusion 
  "Temperature (Excluded)" = "#BBBBBB",
  "Temperature (Included)" = "orange2"
)

# produce the plot showing all data points over time grouped based on inclusion vs exclusion 
temp_plot <- ggplot(dat_ex_temp, aes(x = month, y = temp)) +
  # annotate boundaries of data inclusion
  geom_vline(xintercept = 44, linewidth = 1.5, colour = "#999999", linetype = "dotted") +
  annotate("text", x = 43.8, y = max(dat_ex_temp$temp) * 0.99, label = "August 2022", 
           angle = 90, vjust = -0.8, size = 4, color = "#999999", fontface = "bold") +
  geom_vline(xintercept = 67, linewidth = 1.5, colour = "#999999", linetype = "dotted") +
  annotate("text", x = 67.1, y = max(dat_ex_temp$temp) * 0.99, label = "June 2024", 
           angle = 90, vjust = -0.8, size = 4, color = "#999999", fontface = "bold") +
  # plot all three groups as points connected by lines
  geom_line(data = filter(dat_ex_temp, Group == "temp_Excluded_1"),
            aes(color = Group), size = 0.5) +
  geom_line(data = filter(dat_ex_temp, Group == "temp_Included"),
            aes(color = Group), size = 0.5) +
  geom_line(data = filter(dat_ex_temp, Group == "temp_Excluded_2"),
            aes(color = Group), size = 0.5) +
  geom_point(aes(color = Group), size = 3L) +
  # reference line colors
  scale_color_manual(
    values = group_colors_3,
    guide = "none"
  ) +
  
  geom_point(aes(fill = LegendGroup), shape = 21, size = 4, color = NA, show.legend = TRUE) +
  
  scale_fill_manual(
    values = group_colors_2,
    name = "",
    labels = names(group_colors_2)
  ) +
  # define axis names 
  labs(
    x = "Month",
    y = "Degrees Celcius"
  ) +
  theme_minimal(base_size = 14) +
  # visually format plot
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = c(0.25, 0.90),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

temp_plot






