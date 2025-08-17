library(ggplot2)
library(esquisse)
library(sensitivity)
library(epiR)
library(dplyr)
library(lhs)
library(readxl)


source("cholera_dpm_vax_sa.R")

cholera_dpm_vax_safun <- function(parameters) {
  with(as.list(c(parameters)), {
    out <- tryCatch({
      withCallingHandlers(
        {ode(y = istate, times = tps, func = cholera_dpm_vax_sa, parms = parameters)}, 
        warning = function(w) {
          invokeRestart("muffleWarning")
        })}, 
      error = function(e) {
        return(NULL)
      }
    )
    
    if (is.null(out)) {
      return(c(NA_real_, NA_real_))
    }
    
    out_c <- data.frame(out)
    
    if (!all(tps %in% out_c$time)) {
      return(c(NA_real_, NA_real_))
    }
    
    if (any(is.na(out_c$CInc)) || 
        any(is.infinite(out_c$CInc)) ||
        any(is.na(out_c$CDc)) ||
        any(is.infinite(out_c$CInc))) {
      return(c(NA_real_, NA_real_))
    }
    
    c_inc_tail <- tail(out_c$CInc, 1) # total cases
    d_inc_tail <- tail(out_c$CDc, 1) # total deaths
    
    return(c(c_inc_tail,d_inc_tail))
  })
}

parm.ranges = list(
  mu = c(1/(67.8*365),1/(67.8*365)), # NO TEST # human natural birth/death rate - 1/(68.43*365)
  
  beta_e0 = c(2.25e-6,2.25e-2), # TEST # average environmental contact rate
  beta_h = c(4.25e-8,4.25e-4), # TEST # human contact rate
  
  gamma = c(1/5, 1/0.5), # TEST # incubation period 
  j = c(0.02,0.02), # NO TEST (via literature) # proportion of cases severe, no vax
  f = c(0.05,0.05), # NO TEST (via literature) # proportion of cases moderate, no vax
  g = c(0.18,0.18), # NO TEST (via literature) #proportion of cases mild, no vax
  
  omega_cf = c(1/3464.296,1/3464.296), # NO TEST (via literature) # first year loss of natural immunity, clinical case
  omega_cp = c(1/730,1/1460), # TEST (via literature) # after first year loss of natural immunity, clinical case 
  omega_a = c(1/(9*30), 1/(6*30)), # TEST # loss of natural immunity, asymptomatic case
  
  theta_12 = c(1/3,1/3), # NO TEST # untreated severe case improves to moderate
  theta_23 = c(1/3,1/3), # NO TEST # untreated moderate case improves to mild
  theta_34 = c(1/3,1/3), # NO TEST # untreated mild case improves to prev. clinical asymptomatic 
  tau_c = c(1,1), # NO TEST # prev. clinical asymptomatic case makes full recovery
  tau_a = c(1/2,1), # TEST # non-clinical asymptomatic case makes full recovery
  
  tau_t = c(0.73,0.88), # TEST # recovery of treated severe cases
  
  delta_1 = c(0.5,0.5), # NO TEST (via literature) # cholera death rate of untreated severe cases
  delta_2 = c(0.1,0.25), # TEST # cholera death rate of untreated moderate cases
  
  zeta = c(1e9,1e12), # TEST # severe case human bacterial shedding rate
  sigma_2 = c(0.1,0.1), # NO TEST # relative shedding rate for moderate case
  sigma_3 = c(0.01,0.01), # NO TEST # relative shedding rate for mild case 
  sigma_4 = c(0.0001,0.0001), # NO TEST # relative shedding rate for asymptomatic case
  
  rr_cd = c(0.29,0.29), # NO TEST # reporting rate of severe cholera cases who died from cholera
  rr_t = c(1,1), # NO TEST # reporting rate of severe cholera cases who were treated
  
  v_H = c(5,5), # NO TEST # hyperinfectious bacteria natural death rate/state-transition rate
  b = c(0.1,0.1), # NO TEST #c(0.1,1), # TEST # proportion of hyperinfectious bacteria secreted into the natural environment
  
  n_0 = c(0.33,0.33),  # NO TEST (via literature) # mean natural state bacterial intrinsic proliferation 
  v_L = c(0.33,0.33),  # NO TEST (via literature) # natural bacteria natural death rate
  KL = c(1e7,1e7), # NO TEST (via literature) # concentration of hyperinfectious V. cholerae in (natural environment) water that yields 50% chance of infection
  K = c(1e4,1e10), # TEST #c(1e7,1e12), # TEST # natural environment carrying capacity concentration (and human environment?)
  
  temp_mean= c(21.99,21.99),  # NO TEST (via fitting) # mean temperature (celsius)
  temp_max = c(23.98,23.98), # NO TEST (via fitting) # maximum temperature (celsius)
  rain_max = c(1067245,1067245), # NO TEST (via fitting) # maximum rainfall (mm)
  A_temp = c(1.18,1.18), # NO TEST (via fitting) # amplitude temperature equation
  A_rain = c(319879.08,319879.08), # NO TEST (via fitting) # amplitude rainfall equation
  phi_temp = c(7.16,7.16), # NO TEST (via fitting) # phase shift temperature equation
  phi_rain = c(8.87,8.87), # NO TEST (via fitting) # phase shift rainfall equation
  rain_mean = c(414168.02,414168.02), # NO TEST (via fitting) # mean rainfall (mm)
  
  alpha_n = c(1,1.8),  # TEST # natural state bacterial intrinsic proliferation sensitivity to temperature
  alpha_b = c(1,2), # TEST # environmental contact rate sensitivity to rainfall
  lag_tr = c(5,5), # NO TEST #temperature and rainfall lag
  
  omega_v1 = c(1/(2*365),1/(1*365)), # TEST (via literature) # loss of vaccine-induced immunity 1D
  omega_v2 = c(1/(7*365),1/(5*365)), # TEST (via literature) # loss of vaccine-induced immunity 2D
  
  a = c(1,3), # TEST (via literature) # assortativity term
  l1 = c(0.587,0.587), # NO TEST (via literature) # vaccine efficacy 1D
  l2 = c(0.65,0.65), # NO TEST (via literature) # vaccine efficacy 2D
  q = c(0.45,0.45), # NO TEST (via literature) # vaccine risk reduction of severe+moderate cases
  r = c(0.4,0.4),  # NO TEST (via literature) # vaccine risk reduction of mild cases
  
  v1_cov = c(0,0.98), # TEST
  rel_v2_cov = c(0,0.98), # TEST
  
  hyg_efct = c(0.05,0.05), # NO TEST bc they directly impact other parameters we are testing #c(0.05,0.8), # TEST # hygiene intervention effectiveness
  san_efct = c(0.07,0.07), #c(0.07,0.8), # TEST # sanitation intervention effectiveness
  cwa_efct = C(0.5,0.5) #c(0.5,0.9) # TEST # clean water access intervention effectiveness
)

### Sensitivity Analysis :: LHS -----
  ## Observing CInc and CDc

num_parms <- 53 # 21; get rid of random bacterial shit; clinical immunity; add vaccine immunity
num_samples <- 150 #200 #for diagnostics, 10-20 runs/variable, for final plot, 500 # testing 15 parameters now 

lhs_matrix <- randomLHS(num_samples, num_parms)

parm.samples <- matrix(NA, nrow = num_samples, ncol = num_parms)
colnames(parm.samples) <- names(parm.ranges)

for (i in seq_len(num_parms)) {
  min_val <- parm.ranges[[i]][1]
  max_val <- parm.ranges[[i]][2]
  parm.samples[,i] <- lhs_matrix[, i] * (max_val - min_val) + min_val
}

parm.samples.df <- as.data.frame(parm.samples)
parm.samples.df[["c_inc_tail"]] <- 0.0
parm.samples.df[["d_inc_tail"]] <- 0.0


f_name <- "lhs.sims.test.rds"
if (!file.exists(f_name)) {
  start_time <- Sys.time()
  for (rr in 1:nrow(parm.samples)) {
    pp <- parm.samples[rr,]
    parm.samples.df[rr, "c_inc_tail"] <- cholera_dpm_vax_safun(pp)[1]
    parm.samples.df[rr, "d_inc_tail"] <- cholera_dpm_vax_safun(pp)[2]
  }
  end_time <- Sys.time()
  saveRDS(parm.samples.df, file = f_name)
}

lhs.sims <- readRDS(f_name)



parm.names <- setdiff(names(lhs.sims), c("c_inc_tail", "d_inc_tail"))

# output_cinc <- "c_inc_tail"
# output_dinc <- "d_inc_tail"


# lhs_ranges <- data.frame(parameter=character(),
#                      min_cinc=numeric(), max_cinc=numeric(), range_cinc=numeric(),
#                      min_dinc=numeric(), max_dinc=numeric(), range_dinc=numeric(),
#                      stringsAsFactors=FALSE)

for (p in parm.names) {
  if (length(unique(lhs.sims[[p]])) == 1) {
    lhs.sims <- lhs.sims %>% select(-all_of(p))
  }
}

varied_parm.names <- setdiff(names(lhs.sims), c("c_inc_tail", "d_inc_tail"))
# output_cinc <- "c_inc_tail"
# output_dinc <- "d_inc_tail"

# X_all <- as.matrix(lhs.sims[, varied_parm.names])    # LHS parameter samples
# Y_cinc_all <- lhs.sims[[output_cinc]]
# Y_dinc_all <- lhs.sims[[output_dinc]]

prcc_df_cinc <- lhs.sims[!is.na(lhs.sims$c_inc_tail), c(varied_parm.names, "c_inc_tail")] 
prcc_df_dinc <- lhs.sims[!is.na(lhs.sims$d_inc_tail), c(varied_parm.names, "d_inc_tail")]

# usable_samples <- !is.na(Y_cinc_all)
# 
# X <- as.data.frame(X_all[usable_samples, , drop = FALSE])  # preserve matrix format
# Y_cinc <- Y_cinc_all[usable_samples]
# Y_dinc <- Y_dinc_all[usable_samples]

# prcc_result_cinc <- pcc(X, Y_cinc, rank = TRUE, nboot = 1000)  # nboot gives CIs via bootstrap, optional but recommended
# prcc_result_dinc <- pcc(X, Y_dinc, rank = TRUE, nboot = 1000)
# 
# # Print results
# print(prcc_result_cinc)
# print(prcc_result_dinc)

prcc_result_cinc <- epi.prcc(prcc_df_cinc, sided.test = 2)
prcc_result_dinc <- epi.prcc(prcc_df_dinc, sided.test = 2)

    
# 
# prccs <- as.data.frame(prcc_result_cinc$PRCC) # chnage to PCC
# prccs$parameter <- rownames(prccs)
# ggplot(prccs, aes(x = reorder(parameter, abs(original)), y = original)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   coord_flip() +
#   geom_errorbar(aes(ymin = `min. c.i.`, ymax = `max. c.i.`), width = 0.2) +
#   labs(title = "Tornado Plot: Partial Rank Correlation Coefficients",
#        x = "Parameter", y = "PRCC") +
#   theme_minimal()
# 
# 
# 
# prccs_d <- as.data.frame(prcc_result_dinc$PRCC)
# prccs_d$parameter <- rownames(prccs_d)
# ggplot(prccs_d, aes(x = reorder(parameter, abs(original)), y = original)) +
#   geom_bar(stat = "identity", fill = "red3") +
#   coord_flip() +
#   geom_errorbar(aes(ymin = `min. c.i.`, ymax = `max. c.i.`), width = 0.2) +
#   labs(title = "Tornado Plot: Partial Rank Correlation Coefficients",
#        x = "Parameter", y = "PRCC") +
#   theme_minimal()

# (optional) rename columns for clarity

colnames(prcc_result_cinc)[colnames(prcc_result_cinc) == "var"] <- "Parameter"
colnames(prcc_result_cinc)[colnames(prcc_result_cinc) == "est"] <- "PRCC"

# Sort by absolute PRCC for plotting
prcc_result_cinc <- prcc_result_cinc[order(abs(prcc_result_cinc$PRCC)), ]
prcc_result_cinc$Parameter <- factor(prcc_result_cinc$Parameter, levels = prcc_result_cinc$Parameter)

# Now plot
library(ggplot2)
ggplot(prcc_result_cinc, aes(x = Parameter, y = PRCC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "PRCC Tornado Plot (c_inc_tail)",
    x = "Parameter",
    y = "PRCC"
  ) +
  theme_minimal()



colnames(prcc_result_dinc)[colnames(prcc_result_dinc) == "var"] <- "Parameter"
colnames(prcc_result_dinc)[colnames(prcc_result_dinc) == "est"] <- "PRCC"

# Sort by absolute PRCC for plotting
prcc_result_dinc <- prcc_result_dinc[order(abs(prcc_result_dinc$PRCC)), ]
prcc_result_dinc$Parameter <- factor(prcc_result_dinc$Parameter, levels = prcc_result_dinc$Parameter)

# Now plot
library(ggplot2)
ggplot(prcc_result_dinc, aes(x = Parameter, y = PRCC)) +
  geom_bar(stat = "identity", fill = "red3") +
  coord_flip() +
  labs(
    title = "PRCC Tornado Plot (d_inc_tail)",
    x = "Parameter",
    y = "PRCC"
  ) +
  theme_minimal()



### Sensitivity Analysis :: OAT Sampling -----
## Observing CInc and CDc

baseline <- sapply(parm.ranges, function(x) median(x))
oat_steps <- 5 #2 # run with more steps

oat_results <- data.frame(
  parameter = character(), 
  value = numeric(), 
  c_inc_tail = numeric(),
  d_inc_tail = numeric(),
  stringsAsFactors = FALSE
)

for (p in names(parm.ranges)) {
  vals <- seq(parm.ranges[[p]][1], parm.ranges[[p]][2], length.out = oat_steps)
  for (val in vals) {
    test_parms <- baseline
    test_parms[p] <- val
    
    output <- cholera_dpm_vax_safun(test_parms)
    
    oat_results <- rbind(
      oat_results,
      data.frame(
        parameter = p, 
        value = val, 
        c_inc_tail = output[1], 
        d_inc_tail = output[2]
      )
    )
  }
}

print(oat_results)

oat_ranges <- data.frame(
  parameter = character(),
  min_cinc = numeric(), max_cinc = numeric(), range_cinc = numeric(),
  min_dinc = numeric(), max_dinc = numeric(), range_dinc = numeric(),
  stringsAsFactors = FALSE
)

varied_parms <- unique(oat_results$parameter)[
  sapply(unique(oat_results$parameter), function(p) {
    length(unique(oat_results$value[oat_results$parameter == p])) > 1
  })
]

for (p in varied_parms) {
  all_p <- oat_results[oat_results$parameter == p, ]
  min_cinc <- min(all_p$c_inc_tail, na.rm = TRUE)
  max_cinc <- max(all_p$c_inc_tail, na.rm = TRUE)
  range_cinc <- max_cinc - min_cinc
  
  min_dinc <- min(all_p$d_inc_tail, na.rm = TRUE)
  max_dinc <- max(all_p$d_inc_tail, na.rm = TRUE)
  range_dinc <- max_dinc - min_dinc
  
  oat_ranges <- rbind(
    oat_ranges,
    data.frame(
      parameter = p, min_cinc = min_cinc, max_cinc = max_cinc, range_cinc = range_cinc,
      min_dinc = min_dinc, max_dinc = max_dinc, range_dinc = range_dinc
    )
  )
}


print(oat_ranges)


# plot c_inc ranges 
oat_ranges$parameter <- factor(oat_ranges$parameter, levels = oat_ranges$parameter[order(oat_ranges$range_cinc, decreasing = FALSE)])

ggplot(oat_ranges) +
  aes(x = parameter, y = range_cinc) +
  geom_col(fill = "#1465B3") +
  labs(
    x = " ",
    y = "Range of Cumulative Cholera Cases Given Parameter Variations",
    title = "Sensitivity of Key Parameters Observed By Impact of Parameter Variation on Cumulative Cholera Incidence",
    subtitle = "OAT Sampling Method"
  ) +
  coord_flip() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )



# plot d_inc ranges
oat_ranges$parameter <- factor(oat_ranges$parameter, levels = oat_ranges$parameter[order(oat_ranges$range_dinc, decreasing = FALSE)])

ggplot(oat_ranges) +
  aes(x = parameter, y = range_dinc) +
  geom_col(fill = "red3") +
  labs(
    x = " ",
    y = "Range of Cumulative Cholera Deaths Given Parameter Variations",
    title = "Sensitivity of Key Parameters Observed By Impact of Parameter Variation on Cumulative Cholera Deaths",
    subtitle = "OAT Sampling Method"
  ) +
  coord_flip() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
