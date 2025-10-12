# This code was designed to generate plots in the paper.



################
#### set up ####
################

# load SimEngine + functions
{
  library(SimEngine)
  library(ggplot2)
  library(gridExtra)
  library(truncnorm)
  library(survival)
  source("R/create_data.R", local = T)
  source("R/surv_true.R", local = T)
  source("R/surv_km.R", local = T)
  source("R/surv_two.R", local = T)
  # source("R/se_km.R", local = T)
  # source("R/se_two.R", local = T)
  source("R/boot_ci.R", local = T)
}



##############
#### MAIN ####
##############

# It can be ran both locally and in a cluster.

# Fisrt step, for survival function in placebo group, verify the bias and the coverage of it using a Bootstrap method.

# start time
start_time = Sys.time()
set.seed(1018)



n = 10000
surv_time = list(
  "Exp" = list(surv_type = "Exponential", surv_params = 2e-2), # may be some problems
  "Gom" = list(surv_type = "Gompertz", surv_params = c(0.1, 1e-3))
)

# surv_time = list(surv_type = "Gompertz", surv_params = c(0.1, 1e-3))
# surv_time = list(surv_type = "Exponential", surv_params = 2e-2)



p_two_plc = list()
p_two_vac = list()
p_km_plc = list()
p_km_vac = list()
for (k in 1:2) {
  surv = surv_time[[k]]
  
  dat_phaseOne = create_data(n, surv$surv_type, surv$surv_params, "complex")
  dat_phaseTwo_vac = dat_phaseOne %>%
    dplyr::filter(Z == 1 & treat==1) # use phase two data
  dat_phaseOne_plc = dat_phaseOne[dat_phaseOne$treat == 0, ] # treat = 0 in placebo group
  dat_phaseOne_vac = dat_phaseOne[dat_phaseOne$treat == 1, ] # treat = 1 in vaccine group
  model_two_plc = coxph(Surv(Y, delta) ~ X1 + X2, data = dat_phaseOne_plc) # no s in placebo group
  model_two_vac = coxph(Surv(Y, delta) ~ X1 + X2 + S + I(S == 0), data = dat_phaseTwo_vac, weights = ipw) # s in vaccine group
  
  # choose a specific time
  time_max = round(max(dat_phaseOne$Y))
  true_plc = c()
  true_vac = c()
  est_two_plc = c()
  est_two_vac = c()
  est_km_plc = c()
  est_km_vac = c()
  two_plc_low = c()
  two_plc_up = c()
  two_vac_low = c()
  two_vac_up = c()
  km_plc_low = c()
  km_plc_up = c()
  km_vac_low = c()
  km_vac_up = c()
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
    est_km_plc[i] = surv_km(i, dat_phaseOne_plc, "plc")
    est_km_vac[i] = surv_km(i, dat_phaseTwo_vac, "vac")
    surv_ci_plc = boot_ci(dat_phaseOne_plc, i, "plc")
    surv_ci_vac = boot_ci(dat_phaseTwo_vac, i, "vac")
    two_plc_low[i] = surv_ci_plc$two_low
    two_plc_up[i] = surv_ci_plc$two_up
    two_vac_low[i] = surv_ci_vac$two_low
    two_vac_up[i] = surv_ci_vac$two_up
    km_plc_low[i] = surv_ci_plc$km_low
    km_plc_up[i] = surv_ci_plc$km_up
    km_vac_low[i] = surv_ci_vac$km_low
    km_vac_up[i] = surv_ci_vac$km_up
  }
  result = data.frame(time = time_i, true_plc = true_plc, true_vac = true_vac,
                      est_two_plc = est_two_plc, est_two_vac = est_two_vac, est_km_plc = est_km_plc, est_km_vac = est_km_vac,
                      two_plc_low = two_plc_low, two_plc_up = two_plc_up, two_vac_low = two_vac_low, two_vac_up = two_vac_up,
                      km_plc_low = km_plc_low, km_plc_up = km_plc_up, km_vac_low = km_vac_low, km_vac_up = km_vac_up)
  
  # plot the Two-phase sampling survival functions in the placebo group
  p_two_plc[[k]] = ggplot(result, aes(x = time)) +
    geom_line(aes(y = true_plc, color = "true"), linewidth = 1) +
    geom_line(aes(y = est_two_plc, color = "est"), linewidth = 1) +
    geom_ribbon(aes(ymin = two_plc_low, ymax = two_plc_up), fill = "blue", alpha = 0.2) +
    labs(x = "time", y = "survival probability", color = "Legend") +
    theme_minimal() +
    scale_color_manual(values = c("true" = "red", "est" = "blue")) +
    theme(legend.position.inside = c(0.9, 0.9)) +
    ggtitle(paste0("Two-phase sampling survival curves in the placebo group (", surv$surv_type, ")"))
  # plot the Two-phase sampling survival functions in the vaccine group
  p_two_vac[[k]] = ggplot(result, aes(x = time)) +
    geom_line(aes(y = true_vac, color = "true"), linewidth = 1) +
    geom_line(aes(y = est_two_vac, color = "est"), linewidth = 1) +
    geom_ribbon(aes(ymin = two_vac_low, ymax = two_vac_up), fill = "blue", alpha = 0.2) +
    labs(x = "time", y = "survival probability", color = "Legend") +
    theme_minimal() +
    scale_color_manual(values = c("true" = "red", "est" = "blue")) +
    theme(legend.position.inside = c(0.9, 0.9)) +
    ggtitle(paste0("Two-phase sampling survival curves in the vaccine group (", surv$surv_type, ")"))
  
  # plot the Two-phase sampling survival functions in the placebo group
  p_km_plc[[k]] = ggplot(result, aes(x = time)) +
    geom_line(aes(y = true_plc, color = "true"), linewidth = 1) +
    geom_line(aes(y = est_km_plc, color = "est"), linewidth = 1) +
    geom_ribbon(aes(ymin = km_plc_low, ymax = km_plc_up), fill = "blue", alpha = 0.2) +
    labs(x = "time", y = "survival probability", color = "Legend") +
    theme_minimal() +
    scale_color_manual(values = c("true" = "red", "est" = "blue")) +
    theme(legend.position.inside = c(0.9, 0.9)) +
    ggtitle(paste0("K-M survival curves in the placebo group (", surv$surv_type, ")"))
  # plot the Two-phase sampling survival functions in the vaccine group
  p_km_vac[[k]] = ggplot(result, aes(x = time)) +
    geom_line(aes(y = true_vac, color = "true"), linewidth = 1) +
    geom_line(aes(y = est_km_vac, color = "est"), linewidth = 1) +
    geom_ribbon(aes(ymin = km_vac_low, ymax = km_vac_up), fill = "blue", alpha = 0.2) +
    labs(x = "time", y = "survival probability", color = "Legend") +
    theme_minimal() +
    scale_color_manual(values = c("true" = "red", "est" = "blue")) +
    theme(legend.position.inside = c(0.9, 0.9)) +
    ggtitle(paste0("K-M survival curves in the vaccine group (", surv$surv_type, ")"))
}

p_two = grid.arrange(p_two_plc[[1]], p_two_vac[[1]], p_two_plc[[2]], p_two_vac[[2]], ncol = 2)
p_km = grid.arrange(p_km_plc[[1]], p_km_vac[[1]], p_km_plc[[2]], p_km_vac[[2]], ncol = 2)
ggsave("Plots/plots_two.png", plot = p_two, width = 13.3, height = 8.3, dpi = 300)
ggsave("Plots/plots_km.png", plot = p_km, width = 13.3, height = 8.3, dpi = 300)



# end time
end_time = Sys.time()
execution_time = end_time - start_time
print(execution_time)
