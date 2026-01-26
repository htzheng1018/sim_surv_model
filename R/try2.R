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
  source("R/est_tps.r", local = T)
  source("R/est_flx.r", local = T)
}



set.seed(1018)
n = 8000
lambda = c(0.1, 1e-3)



dat = create_data(n, "Gompertz", lambda, "complex", ind = F)
est_tps(dat = dat, t = 40)
est_flx(dat = dat, t = 40)

