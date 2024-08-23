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
dat_phaseTwo = dat_phaseOne %>%
  dplyr::filter(Z == 1 & treat == 1) # use phase two data
wt_phase = nrow(dat_phaseTwo) / nrow(dat_phaseOne)
model = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_phaseTwo, weights = ipw)
summary(model)



est_surv = function(model, t, data) {
  bh_1 = basehaz(model)
  bh_1f = function(t) {
    index = which.min(abs(bh_1$time-t))
    return(bh_1$hazard[index])
  }
  
  coeffs = model$coefficients
  bh_2 = function(data, coeffs) {
    XS = data[, c("X1", "X2", "S")]
    return(exp(coeffs %*% t(XS)))
  }
  
  result = exp(- bh_1f(t) * wt_phase * bh_2(dat_phaseTwo, coeffs))
  return(mean(result))
}

est_surv(model, 12.70725, dat_phaseTwo)


real_surv = function(model, t) {
  Q_fit = survfit(model)
  index = which.min(abs(Q_fit$time - t))
  Q_fit$surv[index]
}
real_surv(model, 46)


t_series = sort(dat_phaseTwo$Y)
t_series
est = c()
real = c()
for (i in 1: length(t_series)) {
  est[i] = est_surv(model, t_series[i], dat_phaseTwo)
  real[i] = real_surv(model, t_series[i])
}
result = data.frame(time = t_series, est = est, real = real)

ggplot(result, aes(x = time)) +
  geom_line(aes(y = est, color = "Estimated"), linewidth = 1) +
  geom_line(aes(y = real, color = "Real"), linewidth = 1) +
  labs(x = "Time", y = "Survival Probability", color = "Legend") +
  theme_minimal() +
  scale_color_manual(values = c("Estimated" = "blue", "Real" = "red")) +
  ggtitle("Survival Curves")
