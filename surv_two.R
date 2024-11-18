# Two Phase Sampling estimated survival function
surv_two = function(model, t, data) {
  bh_1 = basehaz(model, centered = F)
  index = which.min(abs(bh_1$time - t))
  
  Q_0 = bh_1$hazard[index]
  
  beta = model$coefficients
  X_S = data[, c(names(model$coefficients))]
  unprop = exp(beta %*% t(X_S))
  
  result = exp(- Q_0 * unprop)
  return(mean(result))
}