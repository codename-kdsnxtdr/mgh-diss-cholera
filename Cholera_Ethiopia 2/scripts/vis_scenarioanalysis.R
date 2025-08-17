# vis_scenarioanalysis.R ----
#
# Given the prospective analysis outputs calculated in run_scenarioanalysis.R, 
# this script visualizes the prospective analysis results by producing various plots
# Depending on the plot that is being rendered, monthly incidence over time may be 
# calculated using sim_cholera_inc.R or interventional impact relative to baseline 
# scenario (s#3) may be calculated.
#
###

### load relevant packages
library(dplyr)
library(esquisse)
library(tidyr)
library(ggplot2)


### load scenario run results: calculated outputs of interest and base outputs
sa_results <- readRDS("outputs/scenarioanalysis_results.rds")
sa_results <- sa_results %>% arrange(desc(scenario == 3)) # move baseline scenario, scenario 3, to the first row of result dataframe

s1_raw <- readRDS("outputs/scenarioanalysis_outs1.rds")
s2_raw <- readRDS("outputs/scenarioanalysis_outs2.rds")
s3_raw <- readRDS("outputs/scenarioanalysis_outs3.rds")
s4_raw <- readRDS("outputs/scenarioanalysis_outs4.rds")
s5_raw <- readRDS("outputs/scenarioanalysis_outs5.rds")
s6_raw <- readRDS("outputs/scenarioanalysis_outs6.rds")
s7_raw <- readRDS("outputs/scenarioanalysis_outs7.rds")
s8_raw <- readRDS("outputs/scenarioanalysis_outs8.rds")
s9_raw <- readRDS("outputs/scenarioanalysis_outs9.rds")
s10_raw <- readRDS("outputs/scenarioanalysis_outs10.rds")
s11_raw <- readRDS("outputs/scenarioanalysis_outs11.rds")


### produce plots 

## plot.6A ----
#### proportional difference between outputs of interest for the baseline scenario 3 vs. other scenarios 
                  
prop_diff <- data.frame(scenario = c(), prop_diff.inc = c(), prop_diff.d = c(), prop_diff.sminc = c()) # data frame for holding calculated proportional differences

# for each scenario-specific set of outputs, calculate proportional difference from scenario 3 baseline
for (ii in 1:nrow(sa_results)) {
  scenario <- sa_results$scenario[ii]
  
  c_diff.inc <- (sa_results$incidence25_28[ii] - sa_results$incidence25_28[1]) / sa_results$incidence25_28[1] # cumulative incidence prop. diff
  c_diff.d <- (sa_results$deaths25_28[ii] - sa_results$deaths25_28[1]) / sa_results$deaths25_28[1] # cumulative deaths prop. diff
  c_diff.sminc <- (sa_results$sevmod_incidence25_28[ii] - sa_results$sevmod_incidence25_28[1]) / sa_results$sevmod_incidence25_28[1] # cumulative incidence sev/mod prop. diff
 
  c_diffs <- list(scenario = scenario, prop_diff.inc = c_diff.inc, prop_diff.d = c_diff.d, prop_diff.sminc = c_diff.sminc) # combine calculated differences 
  
  prop_diff <- rbind(prop_diff, c_diffs) # add to results dataframe
  
}

prop_diff$scenario <- as.factor(prop_diff$scenario) # define scenario names as factors for plotting 
# convert calculated differences into a percentage
prop_diff$prop_diff.inc <- prop_diff$prop_diff.inc * 100
prop_diff$prop_diff.d <- prop_diff$prop_diff.d * 100
prop_diff$prop_diff.sminc <- prop_diff$prop_diff.sminc * 100 

# pivot data for appropriate plotting 
prop_diff_long <- tidyr::pivot_longer(
  prop_diff,
  cols = c(prop_diff.inc, prop_diff.d, prop_diff.sminc),
  names_to = "outcome_type",
  values_to = "difference"
)

# add scenario numbers as column to be used as factors for plotting 
prop_diff_long$scenario <- factor(
  prop_diff_long$scenario,
  levels = c("3", "1", "2", "4", "5", "6", "7", "8", "9", "10","11") 
)

#
#
# plot proportional differences of all outcomes grouped by scenario compared to the baseline 
prop_diff_long$label <- paste0(round(prop_diff_long$difference, 1), "%") # add a percent sign to values

ggplot(prop_diff_long, aes(x = scenario, y = difference, fill = outcome_type)) + # scenarios on x, difference on y, fill based on difference type
  geom_col(position = "dodge") +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 1),
    vjust = ifelse(prop_diff_long$difference < 0, 1.2, -0.5), 
    color = "black", size = 3
  ) +
  theme_minimal() +
  labs(
    x = "Scenario Number",
    y = "Percent Difference from Baseline",
    fill = ""
  ) +
  scale_fill_manual(
    values = c("prop_diff.inc" = "#66c2a5", "prop_diff.d" = "#fc8d62", "prop_diff.sminc" = "#8da0cb"),
    labels = c(
      "prop_diff.inc" = "Total Incident Cases",
      "prop_diff.d" = "Deaths",
      "prop_diff.sminc" = "Severe/Moderate Incident Cases"
    )
  ) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  annotate( # scenario 3 visualize 0
    "rect",
    xmin = 0.5, xmax = 1.45,   
    ymin = -0.1, ymax = 0,
    alpha = 0.5, fill = "black"
  ) 


# theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
#     axis.title = element_text(face = "bold"),
#     axis.text = element_text(size = 12)
#   )

# #### PLOT 2:: all scenarios, deaths compared to base line 
# prop_diff_deaths <- subset(prop_diff_long, outcome_type == "prop_diff.d")
# prop_diff_deaths$label <- paste0(round(prop_diff_deaths$difference, 1), "%")
# 
# ggplot(prop_diff_deaths, aes(x = scenario, y = difference, fill = outcome_type)) +
#   annotate(
#     "rect",
#     xmin = 0.5, xmax = 1.45,   # adjust for your baseline highlight as desired
#     ymin = -0.1, ymax = 0,
#     alpha = 0.5, fill = "black"
#   ) +
#   geom_col(position = "dodge") +
#   # Add value labels on top of each bar
#   geom_text(
#     aes(label = label),
#     vjust = ifelse(prop_diff_deaths$difference < 0, 1.2, -0.5), # place label slightly above bar for negative bars
#     color = "black",
#     size = 4
#   ) +
#   theme_minimal() +
#   labs(
#     x = "Scenario Number",
#     y = "Percent Difference from Scenario3 Baseline",
#     fill = ""
#   ) +
#   scale_fill_manual(
#     values = c("prop_diff.d" = "#fc8d62"),
#     labels = c("prop_diff.d" = "Deaths")
#   ) +
#   theme(
#     legend.position = c(0.95, 0.05),
#     legend.justification = c("right", "bottom"),
#     legend.background = element_rect(fill = "white", color = NA)
#   )




# ## PLOT 3:: all scenarios, total incidence compared to base line 
# prop_diff_inc <- subset(prop_diff_long, outcome_type == "prop_diff.inc")
# prop_diff_inc$label <- paste0(round(prop_diff_inc$difference, 1), "%")
# 
# ggplot(prop_diff_inc, aes(x = scenario, y = difference, fill = outcome_type)) +
#   annotate(
#     "rect",
#     xmin = 0.5, xmax = 1.45,   # adjust for your baseline highlight as desired
#     ymin = -0.1, ymax = 0,
#     alpha = 0.5, fill = "black"
#   ) +
#   geom_col(position = "dodge") +
#   geom_text(
#     aes(label = label),
#     vjust = ifelse(prop_diff_inc$difference < 0, 1.2, -0.5), # place label slightly above bar for negative bars
#     color = "black",
#     size = 4
#   ) +
#   theme_minimal() +
#   labs(
#     x = "Scenario Number",
#     y = "Percent Difference from Scenario3 Baseline",
#     fill = ""
#   ) +
#   scale_fill_manual(
#     values = c("prop_diff.inc" = "#66c2a5"),
#     labels = c("prop_diff.inc" = "All Incident Cases")
#   ) +
#   theme(
#     legend.position = c(0.95, 0.05),
#     legend.justification = c("right", "bottom"),
#     legend.background = element_rect(fill = "white", color = NA)
#   )



# ## PLOT 4:: all scenarios, total incidence/total sever/mod compared to base line 
# prop_diff_incs <- subset(prop_diff_long, outcome_type != "prop_diff.d")
# prop_diff_incs$label <- paste0(round(prop_diff_incs$difference, 1), "%")
# 
# ggplot(prop_diff_incs, aes(x = scenario, y = difference, fill = outcome_type)) +
#   annotate(
#     "rect",
#     xmin = 0.5, xmax = 1.45,   # adjust for your baseline highlight as desired
#     ymin = -0.1, ymax = 0,
#     alpha = 0.5, fill = "black"
#   ) +
#   geom_col(position = "dodge") +
#   geom_text(
#     aes(label = label),
#     position = position_dodge(width = 1),
#     vjust = ifelse(prop_diff_incs$difference < 0, 1.2, -0.5), # place label slightly above bar for negative bars
#     color = "black",
#     size = 4
#   ) +
#   theme_minimal() +
#   labs(
#     x = "Scenario Number",
#     y = "Percent Difference from Scenario3 Baseline",
#     fill = ""
#   ) +
#   scale_fill_manual(
#     values = c("prop_diff.inc" = "#66c2a5", "prop_diff.sminc" = "#8da0cb"),
#     labels = c("prop_diff.inc" = "All Incident Cases", "prop_diff.sminc" = "Severe/Moderate Incident Cases")
#   ) +
#   theme(
#     legend.position = c(0.95, 0.05),
#     legend.justification = c("right", "bottom"),
#     legend.background = element_rect(fill = "white", color = NA)
#   )





## plot.6B ----

# calculate monthly reported incidences across all scenarios

source("scripts/sim_cholera_inc.R") # source the script for calculation
dat <- read_excel("data/raw_data/Cholera_Data_Monthly_Mule_1.xlsx") # source case data

inc_sim_list <- list()
num_scenarios <- 11  

for (i in 1:num_scenarios) { # for each scenario run sim_cholera_inc.R and store outputs 
  raw_data_name <- paste0("s", i, "_raw") # reference the naming scheme of store raw scenario outputs 
  raw_data <- get(raw_data_name)
  
  CRInc <- raw_data$CRInc
  CRD <- raw_data$CRD
  
  sim_result <- sim_cholera_inc(
    CRInc = CRInc, CRD = CRD,
    c_data = dat$incidence[6:length(dat$incidence)],
    d_data = dat$death[6:length(dat$death)],
    years_init = 4, time_stop = time_stop
  )
  
  colnames(sim_result) <- paste0(colnames(sim_result), ".s", i) # customtize how outputs are saved
  
  inc_sim_list[[i]] <- sim_result
}

inc_sim_cbind <- do.call(cbind, inc_sim_list) # combine all elements in the list column-wise

#
#
# plot the reported monthly incidence of the NCP plan (s1) vs. the baseline (s3) ----
cinc_sim.v.cinc_dat <- ggplot(inc_sim_cbind) +
  geom_line( # scenario 3
    aes(x = month.s3, y = rcaseinc_sim.s3, 
        colour = factor("Scenario 3 (Baseline)", levels = c("Scenario 3 (Baseline)", "Scenario 1 (NCP Plan)", "Data"))),
    linewidth = 0.75
  ) +
  geom_line( # scenario 1
    aes(x = month.s1, y = rcaseinc_sim.s1, 
    colour = factor("Scenario 1 (NCP Plan)", levels = c("Scenario 3 (Baseline)", "Scenario 1 (NCP Plan)", "Data"))),
    linewidth = 1.5
  ) + 
  # geom_point(
  #   aes(x = month.s1, y = rcaseinc_data.s1, 
  #   colour = factor("Data", levels = c("Baseline", "Scenario#1", "Data"))),
  #   size = 3L
  # ) +
  geom_vline(xintercept = 109, 
             linewidth = 1.5,
             colour = "#999999",
             linetype = "dotted") +  # Subtle grey vertical line
  geom_vline(xintercept = 73, 
             linewidth = 1.5,
             colour = "#999999",
             linetype = "dotted") +  # Subtle grey vertical line
  labs(
    x = "Month",
    y = "Reported Monthly Case Incidence",
    #title = "Comparison of Reported and Modeled Case Incidence",
    colour = "Scenario"
  ) +
  scale_colour_manual(values = c(
    "Scenario 3 (Baseline)" = "black",
    "Scenario 1 (NCP Plan)" = "#66c2a5",
    "Data" = "#CC6677")) +
  theme_minimal(base_size = 14) +
  theme(
    #plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = c(0.2, 0.90),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = 109, y = 6000, label = "January 2028", angle = 90, vjust = -0.8,size = 4, color = "#999999", fontface = "bold") +
  annotate("text", x = 73, y = 6000, label = "January 2025", angle = 90, vjust = -0.8, size = 4, color = "#999999", fontface = "bold")

cinc_sim.v.cinc_dat + coord_cartesian(xlim = c(49, 108), ylim = c(0, 6500))

#
#
# plot the reported monthly incidence of the baseline with 2D PREv (s2) vs. the baseline (s3) 
cinc_sim.v.cinc_dat <- ggplot(inc_sim_cbind) +
  geom_line(
    aes(x = month.s3, y = rcaseinc_sim.s3, 
        colour = factor("Scenario 3 (Baseline)", levels = c("Scenario 3 (Baseline)", "Scenario 2"))),
    linewidth = 0.75
  ) +
  geom_line(
    aes(x = month.s2, y = rcaseinc_sim.s2, 
        colour = factor("Scenario 2", levels = c("Scenario 3 (Baseline)", "Scenario 2"))),
    linewidth = 1.5
  ) + 
  geom_vline(xintercept = 109, 
             linewidth = 1.5,
             colour = "#999999",
             linetype = "dotted") +  # Subtle grey vertical line
  geom_vline(xintercept = 73, 
             linewidth = 1.5,
             colour = "#999999",
             linetype = "dotted") +  # Subtle grey vertical line
  labs(
    x = "Month",
    y = "Reported Monthly Case Incidence",
    #title = "Comparison of Reported and Modeled Case Incidence",
    colour = "Scenario"
  ) +
  scale_colour_manual(values = c(
    "Scenario 3 (Baseline)" = "black",
    "Scenario 2" = "purple3"
   )) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = c(0.2, 0.90),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = 109, y = 5500, label = "January 2028", angle = 90, vjust = -0.8,size = 4, color = "#999999", fontface = "bold") +
  annotate("text", x = 73, y = 5500, label = "January 2025", angle = 90, vjust = -0.8, size = 4, color = "#999999", fontface = "bold")

cinc_sim.v.cinc_dat + coord_cartesian(xlim = c(49, 108), ylim = c(0, 6000))

#
#
# plot the reported monthly incidence of limited OCV, full wash (s10) vs. the NCP plan (s1) ----
cinc_sim.v.cinc_dat <- ggplot(inc_sim_cbind) +
  geom_line(
    aes(x = month.s10, y = rcaseinc_sim.s10, 
        colour = factor("Scenario 10", levels = c("Scenario 1 (NCP Plan)", "Scenario 10"))),
    linewidth = 1
  ) + 
  geom_line(
    aes(x = month.s1, y = rcaseinc_sim.s1, 
        colour = factor("Scenario 1 (NCP Plan)", levels = c("Scenario 1 (NCP Plan)", "Scenario 10"))),
    linewidth = 1
  ) +
  geom_vline(xintercept = 109, 
             linewidth = 1.5,
             colour = "#999999",
             linetype = "dotted") +  # Subtle grey vertical line
  geom_vline(xintercept = 73, 
             linewidth = 1.5,
             colour = "#999999",
             linetype = "dotted") +  # Subtle grey vertical line
  labs(
    x = "Month",
    y = "Reported Monthly Case Incidence",
    #title = "Comparison of Reported and Modeled Case Incidence",
    colour = "Scenario"
  ) +
  scale_colour_manual(values = c(
    "Scenario 1 (NCP Plan)" = "#66c2a5",
    "Scenario 10" = "pink3"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = c(0.2, 0.90),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = 109, y = 5500, label = "January 2028", angle = 90, vjust = -0.8,size = 4, color = "#999999", fontface = "bold") +
  annotate("text", x = 73, y = 5500, label = "January 2025", angle = 90, vjust = -0.8, size = 4, color = "#999999", fontface = "bold")

cinc_sim.v.cinc_dat + coord_cartesian(xlim = c(49, 108), ylim = c(0, 6000))

#
#
# plot the reported monthly incidence of the ALL WASH improvement COMBOS (s4-11) vs. the baseline (s3) ----

levels_vec <- c("Scenario 3 (Baseline)", "Scenario 4", "Scenario 5", "Scenario 6", "Scenario 7", "Scenario 8", "Scenario 9", "Scenario 10", "Scenario 11")

cinc_sim.v.cinc_dat <- ggplot(inc_sim_cbind) +
  geom_line(
    aes(x = month.s3, y = rcaseinc_sim.s3, 
      colour = factor("Scenario 3 (Baseline)", levels = levels_vec)),
    size = 0.75
  ) + 
  geom_line(
    aes(x = month.s4, y = rcaseinc_sim.s4, 
        colour = factor("Scenario 4", levels = levels_vec))
  ) +
  geom_line(
    aes(x = month.s5, y = rcaseinc_sim.s5, 
        colour = factor("Scenario 5", levels = levels_vec))
  ) +
  geom_line(
    aes(x = month.s6, y = rcaseinc_sim.s6, 
        colour = factor("Scenario 6", levels = levels_vec))
  ) +
  geom_line(
    aes(x = month.s7, y = rcaseinc_sim.s7, 
        colour = factor("Scenario 7", levels = levels_vec))
  ) +
  geom_line(
    aes(x = month.s8, y = rcaseinc_sim.s8, 
        colour = factor("Scenario 8", levels = levels_vec))
  ) +
  geom_line(
    aes(x = month.s9, y = rcaseinc_sim.s9, 
        colour = factor("Scenario 9", levels = levels_vec))
  ) +
  geom_line(
    aes(x = month.s10, y = rcaseinc_sim.s10, 
        colour = factor("Scenario 10", levels = levels_vec))
  ) +
  geom_line(
    aes(x = month.s11, y = rcaseinc_sim.s11, 
        colour = factor("Scenario 11", levels = levels_vec))
  ) +
  geom_vline(xintercept = 109, 
             linewidth = 1.5,
             colour = "#999999",
             linetype = "dotted") +  # Subtle grey vertical line
  geom_vline(xintercept = 73, 
             linewidth = 1.5,
             colour = "#999999",
             linetype = "dotted") +  # Subtle grey vertical line
  labs(
    x = "Month",
    y = "Reported Monthly Case Incidence",
    #title = "Comparison of Reported and Modeled Case Incidence",
    colour = "Scenario"
  ) +
  scale_colour_manual(values = c(
    "Scenario 3 (Baseline)" = "black",   # original dark teal baseline
    "Scenario 4" = "yellow3", # brighter muted yellow
    "Scenario 5" = "red3", # brighter muted red
    "Scenario 6" = "#1B98E0", # brighter muted blue
    "Scenario 7" = "#D99252", # brighter muted orange
    "Scenario 8" = "#609059", # brighter muted green
    "Scenario 9" = "purple3", # brighter muted purple
    "Scenario 10" = "brown", # brighter muted brown
    "Scenario 11" = "#294e42"
  )) +
  theme_minimal(base_size = 14) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.width = unit(2, "lines"),
    legend.key.height = unit(1, "lines"),
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = 109, y = 5500, label = "January 2028", angle = 90, vjust = -0.8,size = 4, color = "#999999", fontface = "bold") +
  annotate("text", x = 73, y = 5500, label = "January 2025", angle = 90, vjust = -0.8, size = 4, color = "#999999", fontface = "bold")

cinc_sim.v.cinc_dat + coord_cartesian(xlim = c(50, 120), ylim = c(0, 6000))
