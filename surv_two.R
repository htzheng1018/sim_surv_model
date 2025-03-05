# Two Phase Sampling estimated survival function
surv_two = function(model, t, data) {
  bh_func = basehaz(model, centered = F)
  index = which.min(abs(bh_func$time - t))
  H_0 = bh_func$hazard[index]
  beta = model$coefficients
  X_S = data[, c(names(model$coefficients))]
  unprop = exp(beta %*% t(X_S))

  result = exp(- H_0 * unprop)
  return(mean(result))
}