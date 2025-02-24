


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


n = 1000
surv_time = list(
  "Exp" = list(surv_type = "Exponential", surv_params = 1.5e-3), # may be some problems
  "Gom" = list(surv_type = "Gompertz", surv_params = c(0.2138, 7e-8))
)

surv_time = list(surv_type = "Gompertz", surv_params = c(0.1, 1e-3))
# surv_time = list(surv_type = "Exponential", surv_params = 2e-2)



dat_phaseOne = create_data(n, surv_time$surv_type, surv_time$surv_params, "iid")
print(max(dat_phaseOne$Y))
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
for (i in 1: time_max) {
  true_plc[i] = surv_true(surv_time$surv_type, surv_time$surv_params, i, dat_phaseOne, "plc")
  true_vac[i] = surv_true(surv_time$surv_type, surv_time$surv_params, i, dat_phaseOne_vac, "vac")
  est_two_plc[i] = surv_two(model_two_plc, i, dat_phaseOne_plc, "plc")
  est_two_vac[i] = surv_two(model_two_vac, i, dat_phaseTwo_vac, "vac")
}
result = data.frame(time = (1: time_max), true_plc = true_plc, true_vac= true_vac, est_plc = est_two_plc, est_vac = est_two_vac)



surv_true(surv_time$surv_type, surv_time$surv_params, 44, dat_phaseOne_vac, "vac")
surv_two(model_two_vac, 44, dat_phaseTwo_vac, "vac")


# plot the survival function in placebo group
ggplot(result, aes(x = time)) +
  geom_line(aes(y = true_plc, color = "true_plc"), linewidth = 1) +
  geom_line(aes(y = est_two_plc, color = "est_two_plc"), linewidth = 1) +
  labs(x = "Time", y = "Survival Probability", color = "Legend") +
  theme_minimal() +
  scale_color_manual(values = c("true_plc" = "red", "est_two_plc" = "blue")) +
  theme(legend.position = c(0.9, 0.9)) +
  ggtitle("Survival Curves")

# plot the survival function in vaccine group
ggplot(result, aes(x = time)) +
  geom_line(aes(y = true_vac, color = "true_vac"), linewidth = 1) +
  geom_line(aes(y = est_two_vac, color = "est_two_vac"), linewidth = 1) +
  labs(x = "Time", y = "Survival Probability", color = "Legend") +
  theme_minimal() +
  scale_color_manual(values = c("true_vac" = "red", "est_two_vac" = "blue")) +
  theme(legend.position = c(0.9, 0.9)) +
  ggtitle("Survival Curves")



# end time
end_time = Sys.time()
execution_time = end_time - start_time
print(execution_time)
