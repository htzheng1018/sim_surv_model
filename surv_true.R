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
      unprop = 0.5*data$X1 + 0.7*data$X2 # no s in placebo group
    } else if (type == "vac") {
      unprop = 0.5*data$X1 + 0.7*data$X2 - 2*data$S # have s in vaccine group
    }
    result = mean(Q_0 ^ (exp(unprop)))
  } else if (integral == "math") {
    if (type == "plc") {
      integrand = function(u) {
        return(Q_0 ^ exp(u))
      }
      result1 = integrate(integrand, lower = 0, upper = 0.7)$value
      result2 = integrate(integrand, lower = 0.5, upper = 1.2)$value
      result = (5/7)*(result1 + result2)
    } else if (type == "vac") {
      # P(S | X1, X2, treat)
      # treat = 1 in vaccine group
      # S = 0 (to avoid integrating the Dirac function)
      integrand1 = function(X1, X2) {
        unprop = 0.5*X1 + 0.7*X2
        Q = Q_0 ^ exp(unprop)
        prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
        result = Q * prob_tmp
        return(result)
      }
      result00 = integrate(function(X2) integrand1(X1 = 0, X2), lower = 0, upper = 1)$value
      result01 = integrate(function(X2) integrand1(X1 = 1, X2), lower = 0, upper = 1)$value
      # S > 0
      integrand2 = function(X1, X2, S) {
        unprop = 0.5*X1 + 0.7*X2 - 2*S
        Q = Q_0 ^ exp(unprop)
        prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
        result = Q * (1 - prob_tmp) * dtruncnorm(S, a = 0, b = 1, mean = 0.5, sd = 0.2)
        return(result)
      }
      result10 = integral2(function(X2, S) integrand2(X1 = 0, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
      result11 = integral2(function(X2, S) integrand2(X1 = 1, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
      # expectations
      result = 0.5*(result00 + result10) + 0.5*(result01 + result11)
    }
  }
  
  return(result)
  # return(mean(Q_0))
  # return(mean(exp(unprop)))
}