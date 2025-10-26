# true survival function
surv_true = function(surv_type, surv_params, t, data, type, method, ind = FALSE) {
  if (surv_type == "Exponential") {
    lambda = surv_params
    Q_0 = exp(- lambda * t)
  } else if ( surv_type == "Gompertz") {
    alpha = surv_params[1]
    lambda = surv_params[2]
    Q_0 = exp((lambda/alpha) * (1 - exp(alpha*t)))
  }
  if (ind == TRUE) {
    data$`I(S == 0)TRUE` = ifelse(data$S == 0, 1, 0)
  }
  
  if (method == "sample") {
    if (type == "plc") {
      unprop = exp(0.5*data$X1 + 0.7*data$X2) # no s in placebo group
    } else if (type == "vac") {
      if (ind == TRUE) {
        unprop = exp(0.5*data$X1 + 0.7*data$X2 - 2*data$S - 0.5*data$`I(S == 0)TRUE`) # have s in vaccine group
      } else {
        unprop = exp(0.5*data$X1 + 0.7*data$X2 - 2*data$S) # have s in vaccine group
      }
    } else if (type == "med") {
      if (ind == TRUE) {
        unprop = exp(0.5*data$X1 + 0.7*data$X2 - 2*0 - 0.5*1) # s = 0 in mediation group
      } else {
        unprop = exp(0.5*data$X1 + 0.7*data$X2 - 2*0) # s = 0 in mediation group
      }
    }
    result = mean(Q_0 ^ (unprop))
  } else if (method == "math") {
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
      if (ind == TRUE) {
        # S = 0 (to avoid integrating the Dirac function)
        integrand1 = function(X1, X2) {
          unprop = 0.5*X1 + 0.7*X2 - 0.5*1
          Q = Q_0 ^ exp(unprop)
          prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
          result = Q * prob_tmp
          return(result)
        }
        result00 = integrate(function(X2) integrand1(X1 = 0, X2), lower = 0, upper = 1)$value
        result01 = integrate(function(X2) integrand1(X1 = 1, X2), lower = 0, upper = 1)$value
        # S > 0
        integrand2 = function(X1, X2, S) {
          unprop = 0.5*X1 + 0.7*X2 - 2*S - 0.5*0
          Q = Q_0 ^ exp(unprop)
          # Q = pmax(Q_0 ^ exp(unprop), 1e-2)
          prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
          result = Q * (1 - prob_tmp) * dtruncnorm(S, a = 0, b = 1, mean = 0.5, sd = 0.2)
          return(result)
        }
        result10 = integral2(function(X2, S) integrand2(X1 = 0, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
        result11 = integral2(function(X2, S) integrand2(X1 = 1, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
        # expectations
        result = 0.5*(result00 + result10) + 0.5*(result01 + result11)
      } else {
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
          # Q = pmax(Q_0 ^ exp(unprop), 1e-2)
          prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
          result = Q * (1 - prob_tmp) * dtruncnorm(S, a = 0, b = 1, mean = 0.5, sd = 0.2)
          return(result)
        }
        result10 = integral2(function(X2, S) integrand2(X1 = 0, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
        result11 = integral2(function(X2, S) integrand2(X1 = 1, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
        # expectations
        result = 0.5*(result00 + result10) + 0.5*(result01 + result11)
      }
    } else if (type == "med") {
      if (ind == TRUE) {
        integrand = function(u) {
          return(Q_0 ^ exp(u))
        }
        result1 = integrate(integrand, lower = -0.5, upper = 0.2)$value
        result2 = integrate(integrand, lower = 0, upper = 0.7)$value
        result = (5/7)*(result1 + result2)
      } else {
        integrand = function(u) {
          return(Q_0 ^ exp(u))
        }
        result1 = integrate(integrand, lower = 0, upper = 0.7)$value
        result2 = integrate(integrand, lower = 0.5, upper = 1.2)$value
        result = (5/7)*(result1 + result2)
      }
    }
  }
  return(result)
}



# Two Phase Sampling estimated survival function
surv_two = function(model, t, data, type, ind = FALSE) {
  bh_func = basehaz(model, centered = F)
  index = which.min(abs(bh_func$time - t))
  H_0 = bh_func$hazard[index]
  beta = model$coefficients
  if (ind == TRUE) {
    data$`I(S == 0)TRUE` = ifelse(data$S == 0, 1, 0)
  }
  X_S = data[, c(names(model$coefficients))]
  unprop = exp(beta %*% t(X_S))
  Q = exp(- H_0 * unprop)
  
  if (type == "plc") {
    result = mean(Q)
  } else if (type == "vac") {
    result = sum(Q * data$ipw) / sum(data$ipw)
  } else if (type == "med") {
    X_S_med = X_S
    X_S_med$`S` = 0
    if (ind == TRUE) {
      X_S_med$`I(S == 0)TRUE` = 1
    }
    unprop_med = exp(beta %*% t(X_S_med))
    Q_med = exp(- H_0 * unprop_med)
    result = sum(Q_med * data$ipw) / sum(data$ipw)
  }
  return(result)
}



create_data = function(n, surv_type, surv_params, sample_type, ind = FALSE) {
  # id
  id = seq(1, n)
  
  # treatment
  treat = sample(0:1, n, replace = TRUE, prob = c(0.3, 0.7))
  
  # X1
  X1 = rbinom(n = n, size = 1, prob = 0.5)
  
  # X2
  X2 = runif(n = n, min = 0, max = 1)
  
  # biomarker
  S = rtruncnorm(n = n, a = 0, b = 1, mean = 0.5, sd = 0.2) # truncated normal in (0, 1)
  prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + treat)) # to make edge_prob distributed not too extremely
  val_tmp = rbinom(n, prob = prob_tmp, size = 1)
  S = treat * (1 - val_tmp) * S
  # S = (1 - val_tmp) * S
  
  # survival time
  U = runif(n = n)
  if (ind == "TRUE") {
    V = 0.5*X1 + 0.7*X2 - 2*S - 0.5*I(S == 0)
  } else {
    V = 0.5*X1 + 0.7*X2 - 2*S
  }
  if (surv_type == "Exponential") {
    lambda = surv_params
    t = -log(U) / (lambda * exp(V))
  } else if (surv_type == "Gompertz") {
    alpha = surv_params[1]
    lambda = surv_params[2]
    t = 1/alpha * log(1 - (alpha * log(U)) / (lambda * exp(V)))
  }
  
  # censored time
  U = runif(n = n)
  if (surv_type == "Exponential") {
    lambda = surv_params
    C = - log(U) / (lambda * exp(0.5*X1 + 0.7*X2))
  } else if (surv_type == "Gompertz") {
    alpha = surv_params[1]
    # alpha = 0.01
    lambda = surv_params[2]
    C = 1/alpha * log(1 - (alpha * log(U)) / (lambda * exp(0.5*X1 + 0.7*X2)))
  }
  
  # delta
  delta = ifelse(t <= C, 1, 0)
  
  # observed time
  Y = pmin(t, C)
  
  # two-phase indicator
  if (sample_type == "iid") {
    prob_tmp = 0.4
    Z = treat * rbinom(n = n, size = 1, prob = prob_tmp)
  } else if (sample_type == "complex") {
    t0 = 80 # set the time of interest
    prob_tmp = delta*I(Y <= t0) + (1 - delta*I(Y <= t0)) * (1 / (1 + exp(- 0.5*X1 - 0.7*X2 + 3)))
    Z = treat * rbinom(n = n, size = 1, prob = prob_tmp)
  }
  
  # temporary dataframe
  data = data.frame("id" = id, "treat" = treat, "Y" = Y, "delta" = delta, "S" = S, "X1" = X1, "X2" = X2, "Z" = Z,
                    "C" = C, "t" = t)
  
  # using ipwpoint function to generate inverse probability weights
  # if (sample_type == "iid") {
  #   ipw = ifelse(Z == 1, 1 / prob_tmp, NA)
  # } else if (sample_type == "complex") {
  #   ipw = ipwpoint(
  #     exposure = Z,
  #     family = "binomial",  # The treatment is binary
  #     link = "logit",
  #     denominator = ~ X1 + X2 + treat,
  #     data = data
  #   )$ipw.weights
  # }
  ipw = ifelse(Z == 1, 1 / prob_tmp, NA)
  
  # final data
  data = data %>% 
    dplyr::mutate(ipw = ipw)
  
  return(data)
}








