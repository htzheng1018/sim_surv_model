




est_tps = function(dat, t, boots = 1000) {
  # get data in different treatment groups
  dat_plc = dat[dat$treat == 0, ]
  dat_vac = dat[dat$Z == 1 & dat$treat == 1, ] # vaccine group in a two-phase sampling framework
  
  # survival functions
  model_plc = coxph(Surv(Y, delta) ~ X1 + X2, data = dat_plc)
  model_vac = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_vac, weights = ipw)
  
  # estimaton of NIE, NDE and proportion mediated
  risk_n = function(model, dat, t, treatment) {
    bh = basehaz(model, centered = F)
    index = which.min(abs(bh$time - t))
    H_0 = bh$hazard[index]
    beta = model$coefficients
    X_S = dat[, c(names(model$coefficients))]
    unprop = exp(beta %*% t(X_S))
    Q = exp(- H_0 * unprop)
    
    if (treatment == "plc") {
      Q_n = mean(Q)
    } else if (treatment== "vac") {
      Q_n = sum(Q * dat$ipw) / sum(dat$ipw)
    } else if (treatment == "med") {
      X_S_med = X_S
      X_S_med$`S` = 0
      unprop_med = exp(beta %*% t(X_S_med))
      Q_med = exp(- H_0 * unprop_med)
      Q_n = sum(Q_med * dat$ipw) / sum(dat$ipw)
    }
    result = 1 - Q_n
  }
  
  r_p = risk_n(model_plc, dat_plc, t, "plc")
  r_v = risk_n(model_vac, dat_vac, t, "vac")
  r_m = risk_n(model_vac, dat_vac, t, "med")
  NIE_n = r_v / r_m
  NDE_n = r_m / r_p
  PM_n = 1 - log(NDE_n) / (log(NIE_n) + log(NDE_n))
  
  # standard error of NIE, NDE and proportion mediated
  # Bootstrap
  nn = nrow(dat)
  R = boots
  NIE.boot = c()
  NDE.boot = c()
  PM.boot = c()
  for (i in 1: R) {
    est.boot = c()
    samps = sample(1: nn, size = nn, replace = TRUE)
    data.boot = dat[samps, ]
    dat_plc.boot = data.boot[data.boot$treat == 0, ]
    dat_vac.boot = data.boot %>% dplyr::filter(Z == 1 & treat==1)
    model_plc.boot = coxph(Surv(Y, delta) ~ X1 + X2, data = dat_plc.boot)
    model_vac.boot = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_vac.boot, weights = ipw)
    r_p.boot = risk_n(model_plc.boot, dat_plc.boot, t, "plc")
    r_v.boot = risk_n(model_vac.boot, dat_vac.boot, t, "vac")
    r_m.boot = risk_n(model_vac.boot, dat_vac.boot, t, "med")
    NIE.boot[i] = r_v.boot / r_m.boot
    NDE.boot[i] = r_m.boot / r_p.boot
    PM.boot[i] = 1 - log(NDE.boot[i]) / (log(NIE.boot[i]) + log(NDE.boot[i]))
  }
  
  ci_NIE = quantile(NIE.boot, prob = c(0.025, 0.975))
  ci_NDE = quantile(NDE.boot, prob = c(0.025, 0.975))
  ci_PM = quantile(PM.boot, prob = c(0.025, 0.975))
  NIE_se = sqrt(sum((NIE.boot - mean(NIE.boot)) ^ 2) / (R - 1))
  NDE_se = sqrt(sum((NDE.boot - mean(NDE.boot)) ^ 2) / (R - 1))
  PM_se = sqrt(sum((PM.boot - mean(PM.boot)) ^ 2) / (R - 1))
  
  result = data.frame(
    estimate = c(NIE_n, NDE_n, PM_n),
    se = c(NIE_se, NDE_se, PM_se),
    low = c(ci_NIE[1], ci_NDE[1], ci_PM[1]),
    up = c(ci_NIE[2], ci_NDE[2], ci_PM[2])
  )
  rownames(result) = c("NIE", "NDE", "PM")
  return(result)
}


