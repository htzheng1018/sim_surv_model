# 加载必要的库

{
  library(SimEngine)
  library(pracma)
  source("create_data.R", local = T)
  source("surv_true.R", local = T)
  source("surv_km.R", local = T)
  source("surv_two.R", local = T)
  source("se_km.R", local = T)
  source("se_two.R", local = T)
  source("boot_ci.R", local = T)
} 

# start time
start_time = Sys.time()



set.seed(1018)

# 定义参数
surv_type = "Exponential"
surv_params = 2e-2
t = 42

surv_type = "Gompertz"
surv_params = c(0.1, 1e-3)
t = 44

n = 10000



dat_phaseOne = create_data(n, surv_type, surv_params, "complex")

dat_phaseTwo_vac = dat_phaseOne %>%
  dplyr::filter(Z == 1 & treat==1) # use phase two data
# dat_phaseOne_plc = dat_phaseOne[dat_phaseOne$treat == 0, ] # treat = 0 in placebo group
dat_phaseOne_vac = dat_phaseOne[dat_phaseOne$treat == 1, ] # treat = 1 in vaccine group
# model_two_plc = coxph(Surv(Y, delta) ~ X1 + X2, data = dat_phaseOne_plc) # no s in placebo group
model_two_vac = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_phaseTwo_vac, weights = ipw) # s in vaccine group

if (surv_type == "Exponential") {
  # t_plc = 19
  t_vac = 42
} else if (surv_type == "Gompertz") {
  # t_plc = 36
  t_vac = 44
}

# get the Survival probability at the specific time point
# Q_true_plc = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, t_plc, dat_phaseOne, "plc", "math")
Q_true_vac = surv_true(surv_type, surv_params, t_vac, dat_phaseOne_vac, "vac", "math")
# Q_est_two_plc = surv_two(model_two_plc, t_plc, dat_phaseOne_plc)
Q_est_two_vac = surv_two(model_two_vac, t_vac, dat_phaseTwo_vac, "vac")



print(paste("true Q is:", Q_true_vac))
print(paste("estimated Q is:", Q_est_two_vac))




