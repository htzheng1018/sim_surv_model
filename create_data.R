create_data = function(n, surv_type, surv_params) {
  # id
  id = seq(1, n)
  
  # treatment
  treat = sample(0:1, n, replace = TRUE, prob = c(0.3, 0.7))
  
  # X1
  X1 = rnorm(n = n, mean = 24.3, sd = 8.38)
  
  # X2
  X2 = rnorm(n = n, mean = 266.84, sd = 507.82)
  
  # biomarker
  S = rtruncnorm(n = n, a = 0, b = 1, mean = 0.5, sd = 0.2) # truncated normal in (0, 1)
  prob_tmp = 1 / (1 + exp(- 0.05*X1 + 0.005*X2 + treat)) # to make edge_prob distributed not too extremely
  val_tmp = rbinom(n, prob = prob_tmp, size = 1)
  S = treat * ((1 - val_tmp) * S)
  
  # survival time
  U = runif(n = n)
  if (surv_type == "Exponential") {
    lambda = surv_params
    t = -log(U) / (lambda * exp(0.15*X1 + 0.001*X2 - 5*S))
  } else if (surv_type == "Gompertz") {
    alpha = surv_params[1]
    lambda = surv_params[2]
    t = 1/alpha * log(1 - (alpha * log(U)) / (lambda * exp(0.15*X1 + 0.001*X2 - 5*S + treat)))
  }
  
  # censore time
  U = runif(n = n)
  if (surv_type == "Exponential") {
    lambda = surv_params
    C = - log(U) / (lambda * exp(0.15*X1 + 0.001*X2))
  } else if (surv_type == "Gompertz") {
    alpha = surv_params[1]
    lambda = surv_params[2]
    C = (- 1/lambda * log(U) * exp(0.15*X1 + 0.001*X2)) ^ (alpha)
  }
  
  # delta
  delta = ifelse(t <= C, 1, 0)
  
  # observed time
  Y = pmin(t, C)
  
  # two-phase indicator
  # t0 = 20 # set the time of interest
  # prob_tmp = 1 / (1 + exp(0.15*X1 + 0.001*X2 - 1))
  # prob_tmp = 1 / (1 + exp(0.15*X1 + 0.001*X2))
  # Z = delta*I(Y <= t0) + (1 - delta*I(Y <= t0)) * rbinom(n = n, size = 1, prob = prob_tmp)
  prob_tmp = 1
  Z = treat * rbinom(n = n, size = 1, prob = prob_tmp)
  
  # temporary dataframe
  data = data.frame("id" = id, "treat" = treat, "Y" = Y, "delta" = delta, "S" = S, "X1" = X1, "X2" = X2, "Z" = Z)
  
  # using ipwpoint function to generate inverse probability weights
  # ip_weights = ipwpoint(
  #   exposure = Z,
  #   family = "binomial",  # The treatment is binary
  #   link = "logit",
  #   denominator = ~ X1 + X2,
  #   data = data
  # )$ipw.weights
  # ip_weights = ifelse(Z == 1, 1 / prob_tmp, NA)
  # ip_weights = ip_weights / sum(ip_weights, na.rm = TRUE) * length(ip_weights) # normalize
  ip_weights = 1 # !!!!!
  
  # final data
  data = data %>% 
    dplyr::mutate(ipw = ip_weights)
  
  return(data)
}
