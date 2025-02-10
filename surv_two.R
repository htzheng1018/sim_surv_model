# Two Phase Sampling estimated survival function
surv_two = function(model, t, data, type) {
  bh_1 = basehaz(model, centered = F)
  index = which.min(abs(bh_1$time - t))

  H_0 = bh_1$hazard[index]

  beta = model$coefficients
  # if (type == "plc") {
  #   beta = c(0.15, 0.001)
  # } else if (type == "vac") {
  #   beta = c(0.15, 0.001, -5)
  # }
  X_S = data[, c(names(model$coefficients))]
  unprop = exp(beta %*% t(X_S)) 

  result = exp(- H_0 * unprop)
  return(mean(result))
  # return(mean(exp(- H_0)))
  # return(mean(unprop))
  # return(mean(beta[2]))
}