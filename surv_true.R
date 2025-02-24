# true survival function
surv_true = function(surv_type, surv_params, t, data, type) {
  if (surv_type == "Exponential") {
    lambda = surv_params
    Q_0 = exp(- lambda * t)
  } else if ( surv_type == "Gompertz") {
    alpha = surv_params[1]
    lambda = surv_params[2]
    Q_0 = exp((lambda/alpha) * (1 - exp(alpha*t)))
  }
  
  if (type == "plc") {
    unprop = 0.5 * data$X1 + 0.7 * data$X2 # no s in placebo group
  } else if (type == "vac") {
    unprop = 0.5 * data$X1 + 0.7 * data$X2 - 2 * data$S # have s in vaccine group
  }
  result = Q_0 ^ (exp(unprop))
  return(mean(result))
  # return(mean(Q_0))
  # return(mean(exp(unprop)))
}