


est_vac_tps = function(data, t) {
  model = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = data, weights = ipw)
  bh = basehaz(model, centered = F)
  index = which.min(abs(bh$time - t))
  H_0 = bh$hazard[index]
  beta = model$coefficients
  X_S = data[, c(names(model$coefficients))]
  unprop = exp(beta %*% t(X_S))
  Q = exp(- H_0 * unprop)
  result = sum(Q * data$ipw) / sum(data$ipw)
  return(result)
}


