### vis_lhs_prcc.R
#
# This script is used to visualize  outputs from sensitivity analysis using the PRCC
# method applied to parameters sampled generated via LHS.
#
# Taking in output from the run_lhs_prcc.R script in the form of a data frame holding all LHS samples 
# and their calculated outputs and the PRCC results, produce PRCC tornado plots for incidence and death outputs
# and produce boxplots demonstrating model outputs across parameter quartiles
#
###


### load relevant packages
library(ggplot2)
library(gridExtra)


### source outputs from LHS-PRCC
source("scripts/run_lhs_prcc.R") 


### PRCC Tornado Plot (c_inc_2528 and d_inc_2528) ----

# c_inc_2528 :: cumulative incidence over prospective time 
colnames(prcc_result_cinc)[colnames(prcc_result_cinc) == "var"] <- "Parameter" # rename PRCC output columns
colnames(prcc_result_cinc)[colnames(prcc_result_cinc) == "est"] <- "PRCC" # rename PRCC output columns

prcc_result_cinc <- prcc_result_cinc[order(abs(prcc_result_cinc$PRCC)), ] # put outputs in order of the highest incidence PRCC value
prcc_result_cinc$Parameter <- factor(prcc_result_cinc$Parameter, levels = prcc_result_cinc$Parameter) # parameters names become a factor for plotting

# produce incidence tornado plot 
ggplot(prcc_result_cinc, aes(x = Parameter, y = PRCC)) +
  geom_bar(stat = "identity", fill = "#66c2a6") +
  coord_flip() +
  labs(
    #title = "PRCC Tornado Plot (c_inc_2528)",
    x = "Parameter",
    y = "PRCC"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold")
  )


# d_inc_2528 :: cumulative deaths over prospective time 
colnames(prcc_result_dinc)[colnames(prcc_result_dinc) == "var"] <- "Parameter" # rename PRCC output columns
colnames(prcc_result_dinc)[colnames(prcc_result_dinc) == "est"] <- "PRCC" # rename PRCC output columns

prcc_result_dinc <- prcc_result_dinc[order(abs(prcc_result_dinc$PRCC)), ] # put outputs in order of the highest deaths PRCC value
prcc_result_dinc$Parameter <- factor(prcc_result_dinc$Parameter, levels = prcc_result_dinc$Parameter) # parameters names become a factor for plotting

# produce deaths tornado plot
ggplot(prcc_result_dinc, aes(x = Parameter, y = PRCC)) +
  geom_bar(stat = "identity", fill = "#fc8d62") +
  coord_flip() +
  labs(
    #title = "PRCC Tornado Plot (d_inc_2528)",
    x = "Parameter",
    y = "PRCC"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold")
  )



### LHS Boxplots ----

## create a column for each varied parameter that tracks what quartile it is in across each sample 
labels <- c("Q1", "Q2", "Q3", "Q4")

for (param in varied_parm.names) {
  quartile_col <- paste0(param, "_quartile")
  lhs.sims[[quartile_col]] <- cut(lhs.sims[[param]],
                                  breaks = quantile(lhs.sims[[param]], probs = seq(0, 1, 0.25), na.rm = TRUE),
                                  include.lowest = TRUE,
                                  labels = labels
  )
}


## parameters demonstrating significant correlation to c_inc_2528
# beta_e0
beta_e0_c <- ggplot(lhs.sims, aes(x = beta_e0_quartile, y = c_inc_2528)) +
  geom_boxplot(fill = "#66c2a5") +
  labs(
    x = "beta_e0", y = "Cumulative Case Incidence Over Prospective Timeframe"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold")
  )
# K
K_c <- ggplot(lhs.sims, aes(x = K_quartile, y = c_inc_2528)) +
  geom_boxplot(fill = "#66c2a5") +
  labs(
    x = "K", y = ""
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold")
  )
# zeta
zeta_c <- ggplot(lhs.sims, aes(x = zeta_quartile, y = c_inc_2528)) +
  geom_boxplot(fill = "#66c2a5") +
  labs(
    x = "zeta", y = ""
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold")
  )
# beta_h
beta_h_c <- ggplot(lhs.sims, aes(x = beta_h_quartile, y = c_inc_2528)) +
  geom_boxplot(fill = "#66c2a5") +
  labs(
    x = "beta_h", y = ""
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold")
  )

grid.arrange(beta_e0_c, K_c, zeta_c, beta_h_c, ncol = 4)





## parameters demonstrating significant correlation to d_inc_2528
# beta_e0
beta_e0_d <- ggplot(lhs.sims, aes(x = beta_e0_quartile, y = d_inc_2528)) +
  geom_boxplot(fill = "#fc8d62") +
  labs(
    x = "beta_e0", y = "Cumulative Deaths Over Prospective Timeframe"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold")
  )
# tau_t
tau_t_d <- ggplot(lhs.sims, aes(x = tau_t_quartile, y = d_inc_2528)) +
  geom_boxplot(fill = "#fc8d62") +
  labs(
    x = "tau_t", y = ""
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold")
  )
# delta_2
delta_2_d <- ggplot(lhs.sims, aes(x = delta_2_quartile, y = d_inc_2528)) +
  geom_boxplot(fill = "#fc8d62") +
  labs(
    x = "delta_2", y = ""
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold")
  )
# K
K_d <- ggplot(lhs.sims, aes(x = K_quartile, y = d_inc_2528)) +
  geom_boxplot(fill = "#fc8d62") +
  labs(
    x = "K", y = ""
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold")
  )
# zeta
zeta_d <- ggplot(lhs.sims, aes(x = zeta_quartile, y = d_inc_2528)) +
  geom_boxplot(fill = "#fc8d62") +
  labs(
    x = "zeta", y = ""
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold")
  )
# beta_h
beta_h_d <- ggplot(lhs.sims, aes(x = beta_h_quartile, y = d_inc_2528)) +
  geom_boxplot(fill = "#fc8d62") +
  labs(
    x = "beta_h", y = ""
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold")
  )

grid.arrange(beta_e0_d, tau_t_d, delta_2_d, K_d, zeta_d, beta_h_d, ncol = 6)



### EXTRA: Plausible CFR values given varied values of tau_t
cfr_hist <- ggplot(lhs.sims, aes(x = cfr_2528*100)) +
  geom_histogram(
    bins = 15,                # Change bins as needed
    fill = "lightblue",        
    color = "black",        
    alpha = 0.85            
  ) +
  labs(
    #title = "Distribution of CFR Across Range of tau_t",
    x = "Value",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12)
  )

cfr_hist


