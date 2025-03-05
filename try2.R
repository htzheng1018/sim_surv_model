# 加载必要的库
library(pracma)  # 用于数值积分
set.seed(1018)

# 定义参数
# alpha <- 0.1
# lambda <- 1e-3
# t <- 36
lambda = 2e-2
t = 42

# 定义 Q_0
# Q_0 <- exp((lambda / alpha) * (1 - exp(alpha * t)))
Q_0 = exp(- lambda * t)

# 定义 P(S | X1, X2, treat) 函数
# P_S_given_X1_X2_treat <- function(S, X1, X2, treat) {
#   prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + treat))
#   result = ifelse(treat == 0, ifelse(S == 0, 1, 0), # treat = 0: S = 0 with probability 1
#                   ifelse(S == 0, # treat = 1: S = 0 with probability prob_tmp, otherwise truncated normal
#                          prob_tmp * ifelse(S == 0, 1, 0) + (1 - prob_tmp) * dtruncnorm(0, a = 0, b = 1, mean = 0.5, sd = 0.2),
#                          (1 - prob_tmp) * dtruncnorm(S, a = 0, b = 1, mean = 0.5, sd = 0.2)))
#   return(result)
# }

# 计算期望
# treat = 0
# S = 0
integrand00 = function(X1, X2) {
  unprop = 0.5*X1 + 0.7*X2
  Q = Q_0 ^ exp(unprop)
  result = Q
  return(result)
}
result000 = integrate(function(X2) integrand00(X1 = 0, X2), lower = 0, upper = 1)$value # X1 = 0
result001 = integrate(function(X2) integrand00(X1 = 1, X2), lower = 0, upper = 1)$value # X1 = 1

# treat = 1
# S = 0
integrand10 = function(X1, X2) {
  unprop = 0.5*X1 + 0.7*X2
  Q = Q_0 ^ exp(unprop)
  prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
  result = Q * prob_tmp
  return(result)
}
result100 = integrate(function(X2) integrand10(X1 = 0, X2), lower = 0, upper = 1)$value
result101 = integrate(function(X2) integrand10(X1 = 1, X2), lower = 0, upper = 1)$value

# S ~ (0,1)
integrand11 = function(X1, X2, S) {
  unprop = 0.5*X1 + 0.7*X2 - 2*S
  Q = Q_0 ^ exp(unprop)
  prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + 1))
  result = Q * (1 - prob_tmp) * dtruncnorm(S, a = 0, b = 1, mean = 0.5, sd = 0.2)
  return(result)
}
result110 = integral2(function(X2, S) integrand11(X1 = 0, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
result111 = integral2(function(X2, S) integrand11(X1 = 1, X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q

# expectation
result = 0.5*(result100 + result110) + 0.5*(result101 + result111)






# sample
source("create_data.R", local = T)
sf_sample = function(Q_0) {
  n = 1000
  # treat = sample(0:1, n, replace = TRUE, prob = c(0.3, 0.7))
  # X1 = rbinom(n = n, size = 1, prob = 0.5)
  # X2 = runif(n = n, min = 0, max = 1)
  # S = rtruncnorm(n = n, a = 0, b = 1, mean = 0.5, sd = 0.2) # truncated normal in (0, 1)
  # prob_tmp = 1 / (1 + exp(0.5*X1 + 0.7*X2 + treat)) # to make edge_prob distributed not too extremely
  # val_tmp = rbinom(n, prob = prob_tmp, size = 1)
  # S = treat * (1 - val_tmp) * S
  
  
  mydata = create_data(n, "Gompertz", c(0.1, 1e-3), "iid")
  dat_vac = mydata[mydata$treat == 1, ]
  X1 = dat_vac$X1
  X2 = dat_vac$X2
  S = dat_vac$S
  
  unprop = 0.5*X1 + 0.7*X2 - 2*S # have s in vaccine group
  result = mean(Q_0 ^ (exp(unprop)))
  return(result)
}


print(paste("数学期望值:", result))
print(paste("抽样期望值:", sf_sample(Q_0)))



