library(truncnorm)
library(ipw)
library(survival)
library(devtools)
library(parallel)
library(dplyr)
library(ggplot2)



n = 1000
surv_type = "Gompertz"
surv_params = c(0.2138, 7e-8)



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
  S = (1 - val_tmp) * S
  
  # survival time
  U = runif(n = n)
  if (surv_type == "Exponential") {
    lambda = surv_params
    t = -log(U) / (lambda * exp(0.15*X1 + 0.001*X2 - 5*S + treat))
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
  t0 = 200 # set the time of interest
  prob_tmp = 1 / (1 + exp(0.15*X1 + 0.001*X2 - 1))
  Z = delta*I(Y <= t0) + (1 - delta*I(Y <= t0)) * rbinom(n = n, size = 1, prob = prob_tmp)
  
  # temporary dataframe
  data = data.frame("id" = id, "treat" = treat, "Y" = Y, "delta" = delta, "S" = S, "X1" = X1, "X2" = X2, "Z" = Z)
  
  # using ipwpoint function to generate inverse probability weights
  ip_weights = ipwpoint(
    exposure = Z,
    family = "binomial",  # The treatment is binary
    link = "logit",
    denominator = ~ treat + delta + S + X1 + X2,
    data = data
  )$ipw.weights
  #ip_weights = 1 / prob_tmp # we can also use this directly as the ip_weights
  
  # final data
  data = data %>% 
    dplyr::mutate(ipw = ip_weights)
  
  return(data)
}



dat_phaseOne = create_data(n, surv_type, surv_params)
model1 = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_phaseOne)
dat_phaseTwo = dat_phaseOne %>%
  dplyr::filter(Z == 1 & treat == 1) # use phase two data
wt_phase = nrow(dat_phaseTwo) / nrow(dat_phaseOne)
model2 = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_phaseTwo)



surv_est = function(model, t, wt, data) {
  bh_1 = basehaz(model)
  index = which.min(abs(bh_1$time - t))
  
  lambda_0 = bh_1$hazard[index]
  
  beta = model$coefficients
  X_S = data[, c(names(model$coefficients))]
  proportional = exp(beta %*% t(X_S))
  
  result = exp(- lambda_0 * wt * proportional)
  return(mean(result))
}

surv_real = function(model, t) {
  Q_fit = survfit(model)
  index = which.min(abs(Q_fit$time - t))
  Q_fit$surv[index]
}



# bootstrap
boot_ci = function(data, wt, model_phaseone, model_phasetwo) {
  time_max = round(max(data$Y))
  nn = nrow(data)
  R = 1000
  surv.boot = NULL
  
  
  surv_hat = c()
  surv_true = c()
  for (i in 1: time_max) {
    surv_hat[i] = surv_est(model_phasetwo, i, wt, data)
    surv_true[i] = surv_real(model_phaseone, i)
  }

  for (r in 1: R) {
    boot.samp = sample(1: nn, size = nn, replace = TRUE)
    data.boot = data[boot.samp, ]
    model.boot = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = data.boot)
    
    est = c()
    for (i in 1: time_max) {
      est[i] = surv_est(model.boot, i, wt, data.boot)
    }
    surv.boot = rbind(surv.boot, est)
    print(surv.boot)
  }
  ci.boot = apply(surv.boot, 2, quantile, prob = c(0.025, 0.975))
  result = data.frame(time = (1: time_max),
                      low = ci.boot[1, ], up = ci.boot[2, ],
                      hat = surv_hat, true = surv_true)
  return(result)
}

surv_ci = boot_ci(dat_phaseTwo, wt_phase, model_phaseone = model1, model_phasetwo = model2)

plot(surv_ci$time, surv_ci$hat, type = "l", col="red", lwd=3)
lines(surv_ci$time, surv_ci$true, col="green", lwd=3)
lines(surv_ci$time, surv_ci$low, col="blue", lwd=3)
lines(surv_ci$time, surv_ci$up, col="blue", lwd=3)






