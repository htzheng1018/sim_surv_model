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
  source("R/ci.R", local = T)
}



set.seed(1018)
n = 4000
lambda = c(0.1, 1e-3)



dat_phaseOne = create_data(n, "Gompertz", lambda, "complex", ind = T) # phase one data (original)
dat_phaseTwo_vac = dat_phaseOne %>%
  dplyr::filter(Z == 1 & treat==1) # vaccine group in phase two data
dat_phaseOne_plc = dat_phaseOne[dat_phaseOne$treat == 0, ] # treat = 0 in placebo group
dat_phaseOne_vac = dat_phaseOne[dat_phaseOne$treat == 1, ] # treat = 1 in vaccine group

model_two_plc = coxph(Surv(Y, delta) ~ X1 + X2, data = dat_phaseOne_plc) # no s in placebo group
# model_two_vac = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_phaseTwo_vac, weights = ipw)
model_two_vac = coxph(Surv(Y, delta) ~ X1 + X2 + S + I(S == 0), data = dat_phaseTwo_vac, weights = ipw)

t_vac = 42



Q_est_two_vac = surv_two(model_two_vac, t_vac, dat_phaseTwo_vac, "vac", ind = T)
print(Q_est_two_vac)
Q_true_vac = surv_true("Gompertz", lambda, t_vac, dat_phaseOne_vac, "vac", "math", ind = T)
print(Q_true_vac)
Q_est_km_vac = surv_km(t_vac, dat_phaseTwo_vac, "vac")
print(Q_est_km_vac)

Q_est_two_med = surv_two(model_two_vac, t_vac, dat_phaseTwo_vac, "med", ind = T)
print(Q_est_two_med)
Q_true_med = surv_true("Gompertz", lambda, t_vac, dat_phaseOne_vac, "med", "math", ind = T)
print(Q_true_med)
Q_est_km_med = surv_km(t_vac, dat_phaseTwo_vac, "med")
print(Q_est_km_med)



ci_boot_med = ci(dat_phaseTwo_vac, t_vac, "med", "bootstrap", ind = T)
print(ci_boot_med)


