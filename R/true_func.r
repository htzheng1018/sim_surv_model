


true_func = function(surv_type, surv_params, t, data, method, ind = FALSE) {
  if (ind == TRUE) {
    data$`I(S == 0)TRUE` = ifelse(data$S == 0, 1, 0)
  }
  
  if (surv_type == "Exponential") {
    lambda = surv_params
    Q_0 = exp(- lambda * t)
  } else if ( surv_type == "Gompertz") {
    alpha = surv_params[1]
    lambda = surv_params[2]
    Q_0 = exp((lambda/alpha) * (1 - exp(alpha*t)))
  }
  
  if (method == "sample") {
    if (ind == TRUE) {
      unprop_plc = exp(0.5*data$X1 + 0.7*data$X2)
      unprop_vac = exp(0.5*data$X1 + 0.7*data$X2 - 2*data$S - 0.5*data$`I(S == 0)TRUE`)
      unprop_med = exp(0.5*data$X1 + 0.7*data$X2 - 2*0 - 0.5*1)
    } else {
      unprop_plc = exp(0.5*data$X1 + 0.7*data$X2)
      unprop_vac = exp(0.5*data$X1 + 0.7*data$X2 - 2*data$S)
      unprop_med = exp(0.5*data$X1 + 0.7*data$X2 - 2*0)
    }
    r_p = 1 - mean(Q_0 ^ (unprop_plc))
    r_v = 1 - mean(Q_0 ^ (unprop_vac))
    r_m = 1 - mean(Q_0 ^ (unprop_med))
    
    NIE = r_v / r_m
    NDE = r_m / r_p
    lambda = 1 - log(NDE) / (log(NIE) + log(NDE))
  } else if (method == "math") {
    # placebo group
    integrand = function(u) {
      return(Q_0 ^ exp(u))
    }
    Q_1 = integrate(integrand, lower = 0, upper = 0.7)$value
    Q_2 = integrate(integrand, lower = 0.5, upper = 1.2)$value
    Q_p = (5/7)*(Q_1 + Q_2)
    r_p = 1 - Q_p
    # vaccine group
    if (ind == TRUE) {
      # S = 0 (to avoid integrating the Dirac function)
      integrand1 = function(X1, X2) {
        unprop = 0.5*X1 + 0.7*X2 - 0.5*1
        Q = Q_0 ^ exp(unprop)
        prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
        result = Q * prob_tmp
        return(result)
      }
      Q_00 = integrate(function(X2) integrand1(X1 = 0, X2), lower = 0, upper = 1)$value
      Q_01 = integrate(function(X2) integrand1(X1 = 1, X2), lower = 0, upper = 1)$value
      # S > 0
      integrand2 = function(X1, X2, S) {
        unprop = 0.5*X1 + 0.7*X2 - 2*S - 0.5*0
        Q = Q_0 ^ exp(unprop)
        # Q = pmax(Q_0 ^ exp(unprop), 1e-2)
        prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
        result = Q * (1 - prob_tmp) * dtruncnorm(S, a = 0, b = 1, mean = 0.5, sd = 0.2)
        return(result)
      }
      Q_10 = integral2(function(X2, S) integrand2(X1 = 0, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
      Q_11 = integral2(function(X2, S) integrand2(X1 = 1, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
      # expectations
      Q_v = 0.5*(Q_00 + Q_10) + 0.5*(Q_01 + Q_11)
    } else if (ind == FALSE){
      # S = 0 (to avoid integrating the Dirac function)
      integrand1 = function(X1, X2) {
        unprop = 0.5*X1 + 0.7*X2
        Q = Q_0 ^ exp(unprop)
        prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
        result = Q * prob_tmp
        return(result)
      }
      Q_00 = integrate(function(X2) integrand1(X1 = 0, X2), lower = 0, upper = 1)$value
      Q_01 = integrate(function(X2) integrand1(X1 = 1, X2), lower = 0, upper = 1)$value
      # S > 0
      integrand2 = function(X1, X2, S) {
        unprop = 0.5*X1 + 0.7*X2 - 2*S
        Q = Q_0 ^ exp(unprop)
        # Q = pmax(Q_0 ^ exp(unprop), 1e-2)
        prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
        result = Q * (1 - prob_tmp) * dtruncnorm(S, a = 0, b = 1, mean = 0.5, sd = 0.2)
        return(result)
      }
      Q_10 = integral2(function(X2, S) integrand2(X1 = 0, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
      Q_11 = integral2(function(X2, S) integrand2(X1 = 1, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
      # expectations
      Q_v = 0.5*(Q_00 + Q_10) + 0.5*(Q_01 + Q_11)
    }
    r_v = 1 - Q_v
    # mediation group
    if (ind == TRUE) {
      integrand = function(u) {
        return(Q_0 ^ exp(u))
      }
      Q_1 = integrate(integrand, lower = -0.5, upper = 0.2)$value
      Q_2 = integrate(integrand, lower = 0, upper = 0.7)$value
      Q_m = (5/7)*(Q_1 + Q_2)
    } else if (ind == FALSE) {
      integrand = function(u) {
        return(Q_0 ^ exp(u))
      }
      Q_1 = integrate(integrand, lower = 0, upper = 0.7)$value
      Q_2 = integrate(integrand, lower = 0.5, upper = 1.2)$value
      Q_m = (5/7)*(Q_1 + Q_2)
    }
    r_m = 1 - Q_m
    
    NIE = r_v / r_m
    NDE = r_m / r_p
    lambda = 1 - log(NDE) / (log(NIE) + log(NDE))
  }
  
  result = data.frame(
    true = c(NIE, NDE, lambda)
  )
  rownames(result) = c("NIE", "NDE", "Proportion Mediated")
  return(result)
}