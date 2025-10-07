# load SimEngine + functions
{
  library(SimEngine)
  library(survival)
  library(parallel)
  library(truncnorm)
  library(devtools)
  library(ipw)
  library(pracma)
  source("R/create_data.R", local = T)
  source("R/surv_true.R", local = T)
  source("R/surv_km.R", local = T)
  source("R/surv_two.R", local = T)
  source("R/se_km.R", local = T)
  source("R/se_two.R", local = T)
  source("R/boot_ci.R", local = T)
}



set.seed(1018)
n = 1000
lambda = 2e-2



dat_phaseOne = create_data(n, "Exponential", lambda, "complex") # phase one data (original)

dat_phaseOne_plc = dat_phaseOne[dat_phaseOne$treat == 0, ] # treat = 0 in placebo group
model_two_plc = coxph(Surv(Y, delta) ~ X1 + X2, data = dat_phaseOne_plc) # no s in placebo group

model_two_plc = coxph(Surv(Y, delta) ~ X1 + X2, data = dat_phaseOne_plc)



t_plc = 19



Q_true_plc = surv_true("Exponential", lambda, t_plc, dat_phaseOne, "plc", "math")
print(Q_true_plc)

Q_true_plc = surv_true("Exponential", lambda, t_plc, dat_phaseOne, "plc", "sample")
print(Q_true_plc)

Q_est_two_plc = surv_two(model_two_plc, t_plc, dat_phaseOne_plc, "plc")
print(Q_est_two_plc)
