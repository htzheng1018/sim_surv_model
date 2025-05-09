# load SimEngine + functions
{
  library(SimEngine)
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



dat_phaseOne = create_data(n, "Exponential", lambda, "complex")
dat_phaseTwo_vac = dat_phaseOne %>%
  dplyr::filter(Z == 1 & treat==1) # vaccine group in phase two data
dat_phaseOne_med = dat_phaseOne[dat_phaseOne$treat == 1, ]
dat_phaseOne_med$S = 0
dat_phaseTwo_med = dat_phaseTwo_vac # treat = 1 in intermediate group
dat_phaseTwo_med$S = 0 # s = 0 in intermediate group

model_two_med = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_phaseTwo_med) # s =0 in intermediate group
# print(model_two_med)

Q_true_med = surv_true("Exponential", lambda, 40, dat_phaseOne_med, "med", "math")
print(Q_true_med)

Q_est_two_med = surv_two(model_two_med, 40, dat_phaseTwo_med, "med")
print(Q_est_two_med)
