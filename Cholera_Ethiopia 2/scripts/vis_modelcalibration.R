### vis_modelcalibration.R
#
# This script takes in outputs from the base model (cholera_dpm_vax.R) and produces plots used 
# main text and supplementary figures
#
# Naming scheme matches figure name in write-up, some tabs contain multiple figures (different panels or multi plot figures)
# Base model parameters run most of the plots, if chnages needed to prodice relevant plot, they are specified
#
###


### load relevant packages ----
library(gridExtra)
library(ggplot2)
library(tidyr)
library(dplyr)


### run cholera model and calculate relevant outputs
  # including monthly reported case incidence, deaths and CFR
source("scripts/cholera_dpm_vax.R") 


#### produce figures

## plot.5A ----
#### calibrated model output against national case incidence data (January2023-June2024)
##### demonstrating success of model calibration
#
#
model_calib_cinc <- ggplot(inc_sim.base) +
  geom_line( # plot model output as a line
    aes(x = month, y = rcaseinc_sim, color = "Model Output"),
    linewidth = 0.75
  ) + 
  geom_point( # plot data as points
    aes(x = month, y = rcaseinc_data, color = "Data"),
    size = 3L
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "Model Output" = "black",  
      "Data" = "#379c78")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Month",
    y = "Reported Monthly Case Incidence",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = c(0.2,0.89),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

# plot model output within data time frame
model_calib_cinc + coord_cartesian(xlim = c(48, 67), ylim = c(0, 5750))



## plot.S5 ----
### calibrated model output against national death incidence data (January2023-June2024)
##### paneled with CFR output vs. CFR data
##### demonstrating success of model calibration
#
#
# create plot for data vs.model output, specific to deaths 
model_calib_dinc <- ggplot(inc_sim.base) +
  geom_line( # plot data as a line to improve visibility
    aes(x = month, y = rdeathinc_data, color = "Data"),
    size = 0.5,  linetype = "dashed"
  ) +
  geom_point( # plot data as points
    aes(x = month, y = rdeathinc_data, color = "Data"),
    size = 3L
  ) +
  geom_line( # plot model output as a line
    aes(x = month, y = rdeathinc_sim, color = "Model Output"),
    linewidth = 0.75
  ) + 
  scale_color_manual( # differentiate the plots by color
    values = c(
      "Model Output" = "#371700",    
      "Data" = "#d4561e")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Month",
    y = "Reported Monthly Deaths",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = c(0.3,0.89),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

# plot model output within data time frame
model_calib_dinc_zoom <- model_calib_dinc + coord_cartesian(xlim = c(48, 67), ylim = c(0, 80))

#
#
# create plot for data vs.model output, specific to calculated CFR 
model_calib_cfr <- ggplot(inc_sim.base) +
  geom_line( # plot data as a line to improve visibility
    aes(x = month, y = cfr_data*100, color = "Data"),
    size = 0.5,  linetype = "dashed"
  ) +
  geom_point( # add data points
    aes(x = month, y = cfr_data*100, color = "Data"),
    size = 3L
  ) +
  geom_line( # plot model output as a line
    aes(x = month, y = cfr_sim*100, color = "Model Output"),
    linewidth = 0.75
  ) + 
  scale_color_manual( # differentiate the plots by color
    values = c(
      "Model Output" = "black",    
      "Data" = "lightblue3")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Month",
    y = "Case Fatality Ratio (%)",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    # format legend
    legend.position = c(0.3,0.89),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

# plot model output within data time frame
model_calib_cfr_zoom <- model_calib_cfr + coord_cartesian(xlim = c(48, 67), ylim = c(0, 3))

# combine plots into a single panel
grid.arrange(model_calib_cfr_zoom, model_calib_dinc_zoom, ncol=2)




## plot.S21 ----
### side-by-side plots of national case and death incidence data from Ethiopia (August2022-June2024)
##### highlighting data point exclusion from August 2022 to January 2023
#
#
# format data
dat_long <- pivot_longer( # differentiate two data types: incidence and death
  dat,
  cols = c(incidence, death),
  names_to = "Type",
  values_to = "Count"
) %>%
  mutate( # set criteria for data that was included/excluded
    Status = ifelse(Month <= 5, "Excluded", "Included"), # if data earlier than January 2023, label as excluded
    Group = paste(Type, Status, sep = "_") # define new groups for each data type, given included or excluded
  )

group_colors <- c( # set custom colors for inclusion/exclusion groups for both data types
  "incidence_Excluded" = "#BBBBBB",  
  "incidence_Included" = "#379c78",  
  "death_Excluded"     = "#BBBBBB",  
  "death_Included"     = "#d4561e"
)

incidence_data <- dat_long %>% filter(Type == "incidence") # separate data type factors to plot individually
death_data     <- dat_long %>% filter(Type == "death") # separate data type factors to plot individually

#
#
# create plot for incidence data, highlighting excluded points
incidence_plot <- ggplot(incidence_data, aes(x = Month, y = Count, color = Group)) +
  geom_vline(xintercept = 6, linewidth = 1.5, colour = "#999999", linetype = "dotted") +
  annotate("text", x = 6, y = 4000, label = "January 2023", 
           angle = 90, vjust = -0.8, size = 4, color = "#999999", fontface = "bold") +
  geom_line(size = 0.5) +
  geom_point(size = 3L) +
  scale_color_manual(
    values = group_colors,
    labels = c(
      "incidence_Excluded" = "Incidence (Excluded)",
      "incidence_Included" = "Incidence (Included)"
    ),
    name = ""
  ) +
  labs(
    x = "Month",
    y = "Count"#,
    #title = "Incidence Data (Excluded Points Highlighted)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = c(0.89,0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )


#
#
# create plot for death data, highlighting excluded points 
death_plot <- ggplot(death_data, aes(x = Month, y = Count, color = Group)) +
  geom_vline(xintercept = 6, linewidth = 1.5, colour = "#999999", linetype = "dotted") +
  annotate("text", x = 6, y = 58, label = "January 2023", 
           angle = 90, vjust = -0.8, size = 4, color = "#999999", fontface = "bold") +
  geom_line(size = 0.5) +
  geom_point(size = 3L) +
  scale_color_manual(
    values = group_colors,
    labels = c(
      "death_Excluded" = "Death (Excluded)",
      "death_Included" = "Death (Included)"
    ),
    name = ""
  ) +
  labs(
    x = "Month",
    y = ""#,
    #title = "Death Data (Excluded Points Highlighted)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = c(0.89,0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

# arrange incidence plot and death plot side-by-side
grid.arrange(incidence_plot, death_plot, ncol=2)


## plot.S22 ----
#### Beta_e and eta(n) over time, highlighting where beta_e0 and eta_0 fall on the graph relative to oscillations 
##### demonstrate the impact of seasonality relative to the base parameter value 
##### demonstrate the direct and indirect impact of seasonality on transmission regulators 
#
#
# plot of beta_e over time demonstrating seasonal oscillations (directly from rainfall)
betae <- ggplot(out_c) +
  geom_line(
    aes(x = time, y = beta_e, color = "beta_e"),
    linewidth = 1
  ) +
  geom_hline(
    aes(yintercept = beta_e0, color = "beta_e0"),
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "beta_e" = "lightgreen",  
      "beta_e0" = "black")
  ) +
    labs( # add x and y axis labels and a legend
      x = "Timestep",
      y = "Value",
      color = "Legend"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      # format axis text
      axis.title = element_text(face = "bold"),
      # format legend
      legend.position = "top",
      legend.background = element_rect(fill = "white", color = NA),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      # add gridlines
      panel.grid.major = element_line(colour = "#D3D3D3"),
      panel.grid.minor = element_blank()
    )
#
#
# plot it eta over time demonstrating seasonal oscillations (directly from temperature)
eta <- ggplot(out_c) +
  geom_line(
    aes(x = time, y = n, color = "eta (n)"),
    linewidth = 1
  ) +
  geom_hline(
    aes(yintercept = 0.33, color = "eta_0 (n_0)"),
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "eta (n)" = "orange1",  
      "eta_0 (n_0)" = "black")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Timestep",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

#
#
# plot of lam_e over time to visualize seasonal oscillations (indirect from rainfall)
lam_e <- ggplot(out_c) +
  geom_line(
    aes(x = time, y = lam_e, color = "lam_e"),
    linewidth = 1
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "lam_e" = "darkgreen")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Timestep",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

#
#
# plot of natural-state bacterial population to visualize seadsonal oscillations (indirect from temperature)
BL <- ggplot(out_c) +
  geom_line(
    aes(x = time, y = BL, color = "BL"),
    linewidth = 1
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "BL" = "brown")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Timestep",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

# combine all plots into one 
grid.arrange(betae, lam_e, eta, BL, ncol = 2)




## plot.S41----
## Hypothetical WASH coverage campaigns impact showing coverage over time 
#### hyg_cov, san_cov, cwa_cov improving from baseline to NCP goal over arbitrary amount of time
#### demonstrate functionality of functions governing WASH inprovement campaigns 
#
# for correct plots, turn off forcing terms and implement coverage improvement camapigns for all WASH components
#
WASH_func <- ggplot(out_c) +
  geom_line( # line for hygiene
    aes(x = time, y = hyg_cov, color = "hyg_cov"),
    linewidth = 1
  ) +
  geom_line( # line for sanitaton 
    aes(x = time, y = san_cov, color = "san_cov"),
    linewidth = 1
  ) +
  geom_line( # line for clean water access
    aes(x = time, y = cwa_cov, color = "cwa_cov"),
    linewidth = 1
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "hyg_cov" = "yellow",
      "san_cov" = "red",
      "cwa_cov" = "blue")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Timestep",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

WASH_func


## rho1 and rho2 overtime, given 16 historical OCV campaigns (2019 to 2024)
#### see if rho1 and rho2 are tiggered appropriately 
#### demonstrate functionality of vaccination function
#
# plot for rho1
rho1_func <- ggplot(out_c) +
  geom_line(
    aes(x = time, y = rho1.rho1, color = "rho1"),
    linewidth = 1
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "rho1" = "pink3")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Timestep",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

rho1_func <- rho1_func + coord_cartesian(xlim = c(0, 2000), ylim = c(0, 0.03))

#
#
# plot for rho2
rho2_func <- ggplot(out_c) +
  geom_line(
    aes(x = time, y = rho2.rho2, color = "rho2"),
    linewidth = 1
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "rho2" = "purple3")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Timestep",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

rho2_func <- rho2_func + coord_cartesian(xlim = c(0, 2000), ylim = c(0, 1.2))


# combine into one plot
grid.arrange(rho1_func, rho2_func, ncol = 2)



## plot.S42----
## V1 and V2 compartments over time, given historical vaccination campaigns
#### show the number of individuals in each compartment over time 
#### demonstrate impact of rho rate triggers
#
# 
#
vax_impact <- ggplot(out_c) +
  geom_line(
    aes(x = time, y = V1, color = "V1"),
    linewidth = 1
  ) +
  geom_line(
    aes(x = time, y = V2, color = "V2"),
    linewidth = 1
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "V1" = "pink3",
      "V2" = "purple3")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Timestep",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

vax_impact


## Improvement of WASH coverage impact on transmission parameters
#### implement arbitrary improvement of WASH components to goal coverage and show impacted transmission parameters
#### demonstrate functionality of functions governing WASH inprovement campaigns 
#
# for correct plots, turn off forcing terms and implement coverage improvement campaigns for individual WASH components (start at t = 1000, end = 1001, NCP goal cov)
#
# hygiene imapacting lam_h 
hyg_impact_h <- ggplot(out_c) +
  geom_line(
    aes(x = time, y = lam_h, color = "lam_h"),
    linewidth = 1
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "lam_h" = "yellow3")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Timestep",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

hyg_impact_h <- hyg_impact_h + coord_cartesian(xlim = c(500, 1500))

#
#
# hygiene impacting lam_e
hyg_impact_e <- ggplot(out_c) +
  geom_line(
    aes(x = time, y = lam_e, color = "lam_e"),
    linewidth = 1
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "lam_e" = "orange3")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Timestep",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

hyg_impact_e <- hyg_impact_e + coord_cartesian(xlim = c(500, 1500))

grid.arrange(hyg_impact_h, hyg_impact_e, ncol=1)

#
#
# sanitation impact K, which indirectly impacts BL/KL+BL
san_impact_K <- ggplot(out_c) +
  geom_line(
    aes(x = time, y = BL/(KL+BL), color = "BL/(KL+BL)"),
    linewidth = 1
  ) +
  geom_hline(
    aes(yintercept = median(BL[1000:1500]/(BL[1000:1500]+KL[1000:1500]))),
    linetype = "dashed", linewidth = 1
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "BL/(KL+BL)" = "red3")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Timestep",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = 550, y = 0.53, label = "0.538", angle = 0, vjust = -0.8, size = 4, color = "black", fontface = "bold")

#
#
# sanitation indirectly impact BL, which impacts BL/KL+BL term, which impacts lam_e
san_impact_e <- ggplot(out_c) +
  geom_line(
    aes(x = time, y = lam_e, color = "lam_e"),
    linewidth = 1
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "lam_e" = "red4")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Timestep",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

san_impact_K <- san_impact_K + coord_cartesian(xlim = c(500, 1500))
san_impact_e <- san_impact_e + coord_cartesian(xlim = c(500, 1500))

grid.arrange(san_impact_K, san_impact_e, ncol=1)

#
#
# clean water access impacts lam_e
cwa_impact <- ggplot(out_c) +
  geom_line(
    aes(x = time, y = lam_e, color = "lam_e"),
    linewidth = 1
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "lam_e" = "blue3")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Timestep",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

cwa_impact + coord_cartesian(xlim = c(500, 1500))






## plot.S6----
## Ratio of lam_e to lam_h over time 
#### demonstrate significance of lam_e relative to lam_h
#
# for correct plots, turn off forcing terms and turn of WASH campaigns
#
lam_ratio <- ggplot(out_c) +
  geom_line(
    aes(x = time, y = lam_e/lam_h, color = "lam_e/lam_h"),
    linewidth = 2.5
  ) +
  scale_color_manual( # differentiate the plots by color 
    values = c(
      "lam_e/lam_h" = "green4")
  ) +
  labs( # add x and y axis labels and a legend
    x = "Timestep",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # format axis text
    axis.title = element_text(face = "bold"),
    # format legend
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # add gridlines
    panel.grid.major = element_line(colour = "#D3D3D3"),
    panel.grid.minor = element_blank()
  )

lam_ratio + coord_cartesian(xlim = c(500, 1500))


