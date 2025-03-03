# 加载必要的库
library(pracma)  # 用于数值积分

# 定义参数
alpha <- 0.1
lambda <- 1e-3
t <- 40

# 定义 Q_0
Q_0 <- exp((lambda / alpha) * (1 - exp(alpha * t)))

# 定义 unprop 函数
unprop <- function(X1, X2, S, treat) {
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
                         prob_tmp * ifelse(S == 0, Inf, 0) + (1 - prob_tmp)*dtruncnorm(0, a = 0, b = 1, mean = 0.5, sd = 0.2),
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
# 由于 X1 和 treat 是离散的，我们分别计算 X1 = 0/1 和 treat = 0/1 的情况
expectation_X1_0_treat_0 <- integral2(function(X2, S) integrand(X2, S, X1 = 0, treat = 0), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
expectation_X1_0_treat_1 <- integral2(function(X2, S) integrand(X2, S, X1 = 0, treat = 1), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
expectation_X1_1_treat_0 <- integral2(function(X2, S) integrand(X2, S, X1 = 1, treat = 0), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
expectation_X1_1_treat_1 <- integral2(function(X2, S) integrand(X2, S, X1 = 1, treat = 1), xmin = 0, xmax = 1, ymin = 0, ymax = 1)

# 最终期望
expectation <- 0.5 * (0.3 * expectation_X1_0_treat_0$Q + 0.7 * expectation_X1_0_treat_1$Q) +
  0.5 * (0.3 * expectation_X1_1_treat_0$Q + 0.7 * expectation_X1_1_treat_1$Q)




# sample
source("create_data.R", local = T)
sf_sample = function(Q_0) {
  data = create_data(1000, "Gompertz", c(0.1, 1e-3), "iid")
  unprop = 0.5*data$X1 + 0.7*data$X2 - 2*data$S # have s in vaccine group
  result = mean(Q_0 ^ (exp(unprop)))
  return(result)
}


print(paste("数学期望值:", expectation))
print(paste("抽样期望值:", sf_sample(Q_0)))



