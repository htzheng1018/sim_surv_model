# true survival function
surv_true = function(surv_type, surv_params, t, data, type, integral) {
  if (surv_type == "Exponential") {
    lambda = surv_params
    Q_0 = exp(- lambda * t)
  } else if ( surv_type == "Gompertz") {
    alpha = surv_params[1]
    lambda = surv_params[2]
    Q_0 = exp((lambda/alpha) * (1 - exp(alpha*t)))
  }
  
  if (integral == "sample") {
    if (type == "plc") {
      unprop = 0.5 * data$X1 + 0.7 * data$X2 # no s in placebo group
    } else if (type == "vac") {
      unprop = 0.5 * data$X1 + 0.7 * data$X2 - 2 * data$S # have s in vaccine group
    }
    result = mean(Q_0 ^ (exp(unprop)))
  } else if (integral == "math") {
    if (type == "plc") {
      integrand = function(u) {
        return(exp(Q_0 * exp(u))) # math formula of unproportional part in placebo group (no s)
      }
      result = integrate(integrand, lower = 0, upper = 1.2)$value
    } else if (type == "vac") {
      integrand1 = function(X2, S) {
        unprop = 0.7 * X2 - 2 * S
        return(Q_0 ^ exp(unprop))
      }
      integrand2 = function(X2, S) {
        unprop = 0.5 + 0.7 * X2 - 2 * S
        return(Q_0 ^ exp(unprop))
      }
      result1 = integral2(integrand1, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
      result2 = integral2(integrand2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
      result = 0.5 * result1$Q + 0.5 * result2$Q
    }
  }
  
  return(result)
  # return(mean(Q_0))
  # return(mean(exp(unprop)))
}