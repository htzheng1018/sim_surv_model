library(truncnorm)
library(ipw)
library(survival)
library(devtools)
library(parallel)
library(dplyr)
library(ggplot2)



# This file try to verify the bias of Survival estimators, and plot the result of them.
# It can be ran both locally and in a cluster.

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



# get the phase one data
dat_phaseOne = create_data(n, surv_type, surv_params)
model1 = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_phaseOne)
summary(model1)

# get the phase two data
dat_phaseTwo = dat_phaseOne
dat_phaseTwo = dat_phaseOne %>%
  dplyr::filter(Z == 1 & treat == 1) # use phase two data
wt_phase = nrow(dat_phaseTwo) / nrow(dat_phaseOne)
model2 = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_phaseTwo, weights = ipw)
summary(model2)



# true survival function
surv_true = function(surv_type, surv_params, model, t, data) {
  if (surv_type == "Exponential") {
    lambda = surv_params
    Q_0 = exp(- lambda * t)
  } else if ( surv_type == "Gompertz") {
    alpha = surv_params[1]
    lambda = surv_params[2]
    Q_0 = exp(lambda/alpha * (1 - exp(alpha*t)))
  }
  beta = model$coefficients
  X_S = data[, c(names(model$coefficients))]
  result = mean(Q_0 ^ (exp(beta %*% t(X_S))))
  return((result))
}

# Coxph estimated survival function
surv_cox = function(model, t) {
  Q_fit = survfit(model)
  index = which.min(abs(Q_fit$time - t))
  result = Q_fit$surv[index]
  return(result)
}

# Two Phase Sampling estimated survival function
surv_two = function(model, t, wt, data) {
  bh_1 = basehaz(model, centered = F)
  index = which.min(abs(bh_1$time - t))
  
  Q_0 = bh_1$hazard[index]
  
  beta = model$coefficients
  X_S = data[, c(names(model$coefficients))]
  proportional = exp(beta %*% t(X_S))
  
  result = exp(- Q_0 * wt * proportional)
  return(mean(result))
}



# get the complete survival functions of all of them
time_max = round(max(dat_phaseOne$Y))
true = c()
est_two = c()
est_cox = c()
for (i in 1: time_max) {
  true[i] = surv_true(surv_type, surv_params, model1, i, dat_phaseOne)
  est_two[i] = surv_two(model2, i, wt_phase, dat_phaseTwo)
  est_cox[i] = surv_cox(model2, i)
}
result = data.frame(time = (1: time_max), true = true, est = est_two, est_cox = est_cox)



# plot the survival function
ggplot(result, aes(x = time)) +
  geom_line(aes(y = true, color = "true"), linewidth = 1) +
  geom_line(aes(y = est_two, color = "est_twophase"), linewidth = 1) +
  geom_line(aes(y = est_cox, color = "est_cox"), linewidth = 1) +
  labs(x = "Time", y = "Survival Probability", color = "Legend") +
  theme_minimal() +
  scale_color_manual(values = c("true" = "blue", "est_twophase" = "red", "est_cox" = "green")) +
  theme(legend.position = c(0.9, 0.9)) +
  ggtitle("Survival Curves")

