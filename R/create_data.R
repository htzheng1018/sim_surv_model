create_data = function(n, surv_type, surv_params, sample_type) {
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
  # S = treat * (1 - val_tmp) * S
  S = (1 - val_tmp) * S
  
  # survival time
  U = runif(n = n)
  if (surv_type == "Exponential") {
    lambda = surv_params
    t = -log(U) / (lambda * exp(0.5*X1 + 0.7*X2 - 2*S))
  } else if (surv_type == "Gompertz") {
    alpha = surv_params[1]
    lambda = surv_params[2]
    t = 1/alpha * log(1 - (alpha * log(U)) / (lambda * exp(0.5*X1 + 0.7*X2 - 2*S)))
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
