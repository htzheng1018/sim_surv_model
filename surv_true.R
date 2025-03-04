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
      prob_S_cond = function(S, X1, X2, treat) {
        prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + treat))
        result = ifelse(treat == 0, ifelse(S == 0, 1, 0), # treat = 0: S = 0 with probability 1
                        ifelse(S == 0, # treat = 1: S = 0 with probability prob_tmp, otherwise truncated normal
                               prob_tmp * ifelse(S == 0, 1, 0),
                               (1 - prob_tmp) * dtruncnorm(S, a = 0, b = 1, mean = 0.5, sd = 0.2)))
        return(result)
      }
      # define the integrand
      # treat = 0, S = 0
      integrand00 = function(X1, X2) {
        unprop = 0.5*X1 + 0.7*X2
        Q = Q_0 ^ exp(unprop)
        result = Q
        return(result)
      }
      result000 = integrate(function(X2) integrand00(X1 = 0, X2), lower = 0, upper = 1)$value # X1 = 0
      result001 = integrate(function(X2) integrand00(X1 = 1, X2), lower = 0, upper = 1)$value # X1 = 1
      # treat = 1, S = 0
      integrand10 = function(X1, X2) {
        unprop = 0.5*X1 + 0.7*X2
        Q = Q_0 ^ exp(unprop)
        prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
        result = Q * prob_tmp
        return(result)
      }
      result100 = integrate(function(X2) integrand10(X1 = 0, X2), lower = 0, upper = 1)$value
      result101 = integrate(function(X2) integrand10(X1 = 1, X2), lower = 0, upper = 1)$value
      # treat = 1, S > 0
      integrand11 = function(X1, X2, S) {
        unprop = 0.5*X1 + 0.7*X2 - 2*S
        Q = Q_0 ^ exp(unprop)
        prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
        result = Q * (1 - prob_tmp) * dtruncnorm(S, a = 0, b = 1, mean = 0.5, sd = 0.2)
        return(result)
      }
      result110 = integral2(function(X2, S) integrand11(X1 = 0, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
      result111 = integral2(function(X2, S) integrand11(X1 = 1, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
      # expectations
      result = 0.5*(0.3*result000 + 0.7*(result100 + result110)) + 0.5*(0.3*result001 + 0.7*(result101 + result111))
    }
  }
  
  return(result)
  # return(mean(Q_0))
  # return(mean(exp(unprop)))
}