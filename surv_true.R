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
      unprop = function(X1, X2, S) {
        return(0.5*X1 + 0.7*X2 - 2*S)
      }
      survival_function = function(Q_0, unprop) {
        return(Q_0 ^ exp(unprop))
      }
      # P(S | X1, X2, treat)
      prob_S_cond = function(S, X1, X2, treat) {
        prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + treat))
        result = ifelse(treat == 0, ifelse(S == 0, ifelse(S == 0, Inf, 0), 0), # treat = 0: S = 0 with probability 1
                        ifelse(S == 0, # treat = 1: S = 0 with probability prob_tmp, otherwise truncated normal
                               prob_tmp*ifelse(S == 0, Inf, 0) + (1 - prob_tmp) * dtruncnorm(0, a = 0, b = 1, mean = 0.5, sd = 0.2),
                               (1 - prob_tmp) * dtruncnorm(S, a = 0, b = 1, mean = 0.5, sd = 0.2)))
        return(result)
      }
      # define the integrand
      integrand = function(X1, X2, S, treat) {
        unprop = 0.5*X1 + 0.7*X2 - 2*S
        Q = Q_0 ^ exp(unprop)
        return(Q * prob_S_cond(S, X1, X2, treat))
      }
      # expectations
      result00 = integral2(function(X2, S) integrand(X1 = 0, X2, S, treat = 0), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
      result01 = integral2(function(X2, S) integrand(X1 = 0, X2, S, treat = 1), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
      result10 = integral2(function(X2, S) integrand(X1 = 1, X2, S, treat = 0), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
      result11 = integral2(function(X2, S) integrand(X1 = 1, X2, S, treat = 1), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
      result = 0.5*(0.3*result00$Q + 0.7*result01$Q) + 0.5*(0.3*result10$Q + 0.7*result11$Q)
    }
  }
  
  return(result)
  # return(mean(Q_0))
  # return(mean(exp(unprop)))
}