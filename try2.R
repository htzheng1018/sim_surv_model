# 加载必要的库
library(pracma)  # 用于数值积分

# 定义参数
alpha <- 0.1
lambda <- 1e-3
t <- 40

# 定义 Q_0
Q_0 <- exp((lambda / alpha) * (1 - exp(alpha * t)))

# 定义 unprop 函数
unprop <- function(X1, X2, S) {
  return(0.5 * X1 + 0.7 * X2 - 2 * S)
}

# 定义 survival function
survival_function <- function(Q_0, unprop) {
  return(Q_0 ^ exp(unprop))
}

# 定义 P(S | X1, X2, treat) 函数
P_S_given_X1_X2_treat <- function(S, X1, X2, treat) {
  prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + treat))
  result = ifelse(treat == 0, ifelse(S == 0, Inf, 0), # treat = 0: S = 0 with probability 1
                  ifelse(S == 0, # treat = 1: S = 0 with probability prob_tmp, otherwise truncated normal
                         prob_tmp * ifelse(S == 0, Inf, 0) + (1 - prob_tmp) * dtruncnorm(0, a = 0, b = 1, mean = 0.5, sd = 0.2),
                         (1 - prob_tmp) * dtruncnorm(S, a = 0, b = 1, mean = 0.5, sd = 0.2)))
  return(result)
}

# 定义被积函数
integrand <- function(X2, S, X1, treat) {
  unprop_value <- unprop(X1, X2, S, treat)
  survival_value <- survival_function(Q_0, unprop_value)
  return(survival_value * P_S_given_X1_X2_treat(S, X1, X2, treat))
}

# 计算期望
# treat = 0
# S = 0
integrand0 <- function(X1, X2) {
  unprop_value <- unprop(X1, X2, S = 0)
  survival_value <- survival_function(Q_0, unprop_value)
  return(survival_value)
}
expectation_X1_0_treat_0 <- integrate(function(X2) integrand0(X1 = 0, X2), lower = 0, upper = 1)$value
expectation_X1_1_treat_0 <- integrate(function(X2) integrand0(X1 = 1, X2), lower = 0, upper = 1)$value

# treat = 1
# S = 0
integrand1 <- function(X1, X2) {
  unprop_value <- unprop(X1, X2, S = 0)
  survival_value <- survival_function(Q_0, unprop_value)
  prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
  result = survival_value * prob_tmp
  return(result)
}
expectation_X1_0_treat_1_S0 <- integrate(function(X2) integrand1(X1 = 0, X2), lower = 0, upper = 1)$value
expectation_X1_1_treat_1_S0 <- integrate(function(X2) integrand1(X1 = 1, X2), lower = 0, upper = 1)$value

# S ~ (0,1)
integrand2 <- function(X1, X2, S) {
  unprop_value <- unprop(X1, X2, S)
  survival_value <- survival_function(Q_0, unprop_value)
  prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
  result = survival_value * (1 - prob_tmp) * dtruncnorm(S, a = 0, b = 1, mean = 0.5, sd = 0.2)
  return(result)
}
expectation_X1_0_treat_1_S_positive <- integral2(function(X2, S) integrand2(X1 = 0, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
expectation_X1_1_treat_1_S_positive <- integral2(function(X2, S) integrand2(X1 = 1, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q

# expectation
expectation <- 0.5 * (0.3 * expectation_X1_0_treat_0 + 0.7 * (expectation_X1_0_treat_1_S0 + expectation_X1_0_treat_1_S_positive)) +
  0.5 * (0.3 * expectation_X1_1_treat_0 + 0.7 * (expectation_X1_1_treat_1_S0 + expectation_X1_1_treat_1_S_positive))






# sample
source("create_data.R", local = T)
sf_sample = function(Q_0) {
  n = 1000
  treat = sample(0:1, n, replace = TRUE, prob = c(0.3, 0.7))
  X1 = rbinom(n = n, size = 1, prob = 0.5)
  X2 = runif(n = n, min = 0, max = 1)
  S = rtruncnorm(n = n, a = 0, b = 1, mean = 0.5, sd = 0.2) # truncated normal in (0, 1)
  prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + treat)) # to make edge_prob distributed not too extremely
  val_tmp = rbinom(n, prob = prob_tmp, size = 1)
  S = treat * (1 - val_tmp) * S
  
  unprop = 0.5*X1 + 0.7*X2 - 2*S # have s in vaccine group
  result = mean(Q_0 ^ (exp(unprop)))
  return(result)
}


print(paste("数学期望值:", expectation))
print(paste("抽样期望值:", sf_sample(Q_0)))



