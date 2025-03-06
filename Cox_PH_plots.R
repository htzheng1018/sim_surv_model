


################
#### set up ####
################

# load SimEngine + functions
{
  library(SimEngine)
  library(ggplot2)
  source("create_data.R", local = T)
  source("surv_true.R", local = T)
  source("surv_km.R", local = T)
  source("surv_two.R", local = T)
  source("se_km.R", local = T)
  source("se_two.R", local = T)
  source("boot_ci.R", local = T)
}


##############
#### MAIN ####
##############

# It can be ran both locally and in a cluster.

# Fisrt step, for survival function in placebo group, verify the bias and the coverage of it using a Bootstrap method.

# start time
start_time = Sys.time()



n = 8000
surv_time = list(
  "Exp" = list(surv_type = "Exponential", surv_params = 2e-2), # may be some problems
  "Gom" = list(surv_type = "Gompertz", surv_params = c(0.1, 1e-3))
)

# surv_time = list(surv_type = "Gompertz", surv_params = c(0.1, 1e-3))
# surv_time = list(surv_type = "Exponential", surv_params = 2e-2)



p_plc = list()
p_vac = list()
for (k in 1:2) {
  surv = surv_time[[k]]
  
  dat_phaseOne = create_data(n, surv$surv_type, surv$surv_params, "complex")
  dat_phaseTwo_vac = dat_phaseOne %>%
    dplyr::filter(Z == 1 & treat==1) # use phase two data
  dat_phaseOne_plc = dat_phaseOne[dat_phaseOne$treat == 0, ] # treat = 0 in placebo group
  dat_phaseOne_vac = dat_phaseOne[dat_phaseOne$treat == 1, ] # treat = 1 in vaccine group
  model_two_plc = coxph(Surv(Y, delta) ~ X1 + X2, data = dat_phaseOne_plc) # no s in placebo group
  model_two_vac = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_phaseTwo_vac, weights = ipw) # s in vaccine group
  
  # choose a specific time
  time_max = round(max(dat_phaseOne$Y))
  true_plc = c()
  true_vac = c()
  est_two_plc = c()
  est_two_vac = c()
  two_plc_low = c()
  two_plc_up = c()
  two_vac_low = c()
  two_vac_up = c()
  time_i = c()
  for (i in seq(1, time_max, by = 1)) {
    time_i[i] = i
    true_plc[i] = surv_true(surv$surv_type, surv$surv_params, i, dat_phaseOne, "plc", "math")
    tryCatch({
      true_vac[i] = surv_true(surv$surv_type, surv$surv_params, i, dat_phaseOne_vac, "vac", "math")
    }, error = function(e) {
      true_vac[i] = true_vac[i - 1]
      print(i)
      print("error")
    })
    est_two_plc[i] = surv_two(model_two_plc, i, dat_phaseOne_plc, "plc")
    est_two_vac[i] = surv_two(model_two_vac, i, dat_phaseTwo_vac, "vac")
    surv_ci_plc = boot_ci(dat_phaseOne_plc, i, "plc")
    surv_ci_vac = boot_ci(dat_phaseTwo_vac, i, "vac")
    two_plc_low[i] = surv_ci_plc$two_low
    two_plc_up[i] = surv_ci_plc$two_up
    two_vac_low[i] = surv_ci_vac$two_low
    two_vac_up[i] = surv_ci_vac$two_up
  }
  result = data.frame(time = time_i, true_plc = true_plc, true_vac = true_vac, est_plc = est_two_plc, est_vac = est_two_vac,
                      plc_low = two_plc_low, plc_up = two_plc_up, vac_low = two_vac_low, vac_up = two_vac_up)
  
  # plot the survival function in placebo group
  p_plc[[k]] = ggplot(result, aes(x = time)) +
    geom_line(aes(y = true_plc, color = "true_plc"), linewidth = 1) +
    geom_line(aes(y = est_plc, color = "est_plc"), linewidth = 1) +
    geom_ribbon(aes(ymin = plc_low, ymax = plc_up), fill = "blue", alpha = 0.2) +
    labs(x = "Time", y = "Survival Probability", color = "Legend") +
    theme_minimal() +
    scale_color_manual(values = c("true_plc" = "red", "est_plc" = "blue")) +
    theme(legend.position = c(0.9, 0.9)) +
    ggtitle("Survival Curves in the placebo group")
  
  # plot the survival function in vaccine group
  p_vac[[k]] = ggplot(result, aes(x = time)) +
    geom_line(aes(y = true_vac, color = "true_vac"), linewidth = 1) +
    geom_line(aes(y = est_vac, color = "est_vac"), linewidth = 1) +
    geom_ribbon(aes(ymin = vac_low, ymax = vac_up), fill = "blue", alpha = 0.2) +
    labs(x = "Time", y = "Survival Probability", color = "Legend") +
    theme_minimal() +
    scale_color_manual(values = c("true_vac" = "red", "est_vac" = "blue")) +
    theme(legend.position = c(0.9, 0.9)) +
    ggtitle("Survival Curves in the vaccine group")
  
  
}



# end time
end_time = Sys.time()
execution_time = end_time - start_time
print(execution_time)
