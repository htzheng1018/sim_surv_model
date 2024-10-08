library(truncnorm)
library(ipw)
library(survival)
library(devtools)
library(parallel)
library(dplyr)
library(ggplot2)



# This file aims at calculating the standard error of the two-sample estimator, and plot the curves to compare the estimator and the true survival function.
# It can be ran both locally and in a cluster.

n = 1000
surv_type = "Gompertz"
surv_params = c(0.2138, 7e-8)
# surv_type = "Exponential"
# surv_params = 1.5e-2



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
  # prob_tmp = 1 / (1 + exp(0.15*X1 + 0.001*X2 - 1))
  prob_tmp = 0.7
  Z = delta*I(Y <= t0) + (1 - delta*I(Y <= t0)) * rbinom(n = n, size = 1, prob = prob_tmp)
  
  # temporary dataframe
  data = data.frame("id" = id, "treat" = treat, "Y" = Y, "delta" = delta, "S" = S, "X1" = X1, "X2" = X2, "Z" = Z)
  
  # # using ipwpoint function to generate inverse probability weights
  # ip_weights = ipwpoint(
  #   exposure = Z,
  #   family = "binomial",  # The treatment is binary
  #   link = "logit",
  #   denominator = ~ treat + delta + S + X1 + X2,
  #   data = data
  # )$ipw.weights
  ip_weights = ifelse(Z == 1, 1 / prob_tmp, NA)
  ip_weights = ip_weights / sum(ip_weights, na.rm = TRUE) * length(ip_weights) # normalize
  
  # final data
  data = data %>% 
    dplyr::mutate(ipw = ip_weights)
  
  return(data)
}



dat_phaseOne = create_data(n, surv_type, surv_params)
model1 = coxph(Surv(Y, delta) ~ X1 + X2 + S + treat, data = dat_phaseOne)
dat_phaseTwo = dat_phaseOne %>%
  dplyr::filter(Z == 1) # use phase two data
model2 = coxph(Surv(Y, delta) ~ X1 + X2 + S + treat, data = dat_phaseTwo, weights = ipw)
wt_phase = nrow(dat_phaseTwo) / nrow(dat_phaseOne)



# true survival function
surv_true = function(surv_type, surv_params, t, data) {
  if (surv_type == "Exponential") {
    lambda = surv_params
    Q_0 = exp(- lambda * t)
  } else if ( surv_type == "Gompertz") {
    alpha = surv_params[1]
    lambda = surv_params[2]
    Q_0 = exp(lambda/alpha * (1 - exp(alpha*t)))
  }
  unprop = 0.15*data$X1 + 0.001*data$X2 - 5*data$S + data$treat
  result = Q_0 ^ (exp(unprop))
  return(mean(result))
}

# Kaplan Meier estimated survival function
surv_km = function(t, data) {
  km.est = survfit(Surv(Y, delta) ~ 1, data = data, conf.type = "log-log")
  index = which.min(abs(km.est$time - t))
  result = km.est$surv[index]
  return(result)
}

# Two Phase Sampling estimated survival function
surv_two = function(model, t, wt, data) {
  bh_1 = basehaz(model, centered = F)
  index = which.min(abs(bh_1$time - t))
  
  Q_0 = bh_1$hazard[index]
  
  beta = model$coefficients
  X_S = data[, c(names(model$coefficients))]
  unprop = exp(beta %*% t(X_S))
  
  result = exp(- Q_0 * wt * unprop)
  return(mean(result))
}



# get true and estimated survival functions
time_max = round(max(dat_phaseOne$Y))
true = c()
est = c()
for (i in 1: time_max) {
  true[i] = surv_true(surv_type, surv_params, i, dat_phaseOne)
  est[i] = surv_two(model2, i, wt_phase, dat_phaseTwo)
}
result = data.frame(time = (1: time_max), true = true, est = est)

# bootstrap
boot_ci = function(wt, data, time_max) {
  nn = nrow(data)
  R = 100
  surv.boot = NULL
  
  for (r in 1: R) {
    boot.samp = sample(1: nn, size = nn, replace = TRUE)
    data.boot = data[boot.samp, ]
    model.boot = coxph(Surv(Y, delta) ~ X1 + X2 + S + treat, data = data.boot, weights = ipw)
    
    est = c()
    for (i in 1: time_max) {
      est[i] = surv_two(model.boot, i, wt, data.boot)
    }
    surv.boot = rbind(surv.boot, est)
  }
  ci = apply(surv.boot, 2, quantile, prob = c(0.025, 0.975))
  se = sqrt(colSums((surv.boot - colMeans(surv.boot)) ^ 2) / (R - 1))
  result = data.frame(time = (1: time_max),
                      low = ci[1, ], up = ci[2, ], se = se)
  return(result)
}

surv_ci = boot_ci(wt_phase, dat_phaseTwo, time_max)
surv_ci = merge(surv_ci, result, by = "time")

plot(surv_ci$time, surv_ci$est, type = "l", col="red", lwd=3)
lines(surv_ci$time, surv_ci$low, col="green", lwd=3)
lines(surv_ci$time, surv_ci$up, col="green", lwd=3)
lines(surv_ci$time, surv_ci$true, col="blue", lwd=3)





