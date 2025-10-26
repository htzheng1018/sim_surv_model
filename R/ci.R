# Bootstrap function (true survival function and estimators)
ci = function(data, t, type, method, ind = FALSE) {
  if (method == "bootstrap") {
    nn = nrow(data)
    R = 100
    surv_km.boot = c()
    surv_two.boot = c()
    
    if (type == "plc") {
      for (r in 1: R) {
        boot.samp = sample(1: nn, size = nn, replace = TRUE)
        data.boot = data[boot.samp, ]
        model.boot = coxph(Surv(Y, delta) ~ X1 + X2, data = data.boot)
        surv_km.boot[r] = surv_km(t, data.boot, type)
        surv_two.boot[r] = surv_two(model.boot, t, data, type, ind)
      }
    } else if (type == "vac") {
      for (r in 1: R) {
        boot.samp = sample(1: nn, size = nn, replace = TRUE)
        data.boot = data[boot.samp, ]
        if (ind == TRUE) {
          model.boot = coxph(Surv(Y, delta) ~ X1 + X2 + S + I(S == 0), data = data.boot, weights = ipw)
        } else {
          model.boot = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = data.boot, weights = ipw)
        }
        surv_km.boot[r] = surv_km(t, data.boot, type)
        surv_two.boot[r] = surv_two(model.boot, t, data.boot, type, ind)
      }
    } else if (type == "med") {
      for (r in 1: R) {
        boot.samp = sample(1: nn, size = nn, replace = TRUE)
        data.boot = data[boot.samp, ]
        if (ind == TRUE) {
          model.boot = coxph(Surv(Y, delta) ~ X1 + X2 + S + I(S == 0), data = data.boot, weights = ipw)
        } else {
          model.boot = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = data.boot, weights = ipw)
        }
        surv_km.boot[r] = surv_km(t, data.boot, type)
        surv_two.boot[r] = surv_two(model.boot, t, data.boot, type, ind)
      }
    }
  } else if (method == "math") {
    # if (type == "plc") {
    #   model = coxph(Surv(Y, delta) ~ X1 + X2, data = data)
    #   covmat = vcov(model) # covariance matrix for beta and H_0
    #   
    # }
  }
  

  ci_km = quantile(surv_km.boot, prob = c(0.025, 0.975))
  ci_two = quantile(surv_two.boot, prob = c(0.025, 0.975))
  km_se = sqrt(sum((surv_km.boot - mean(surv_km.boot)) ^ 2) / (R - 1))
  two_se = sqrt(sum((surv_two.boot - mean(surv_two.boot)) ^ 2) / (R - 1))
  result = data.frame(km_low = ci_km[1], km_up = ci_km[2], km_se = km_se,
                      two_low = ci_two[1], two_up = ci_two[2], two_se = two_se)
  return(result)
}