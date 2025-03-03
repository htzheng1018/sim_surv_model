# Load necessary library
library(pracma)

# Define parameters
alpha <- 0.1
lambda <- 1e-3
t <- 40

# Define Q_0
Q_0 <- exp((lambda / alpha) * (1 - exp(alpha * t)))

# Define integrand 1
integrand1 <- function(X2, S) {
  unprop <- 0.7 * X2 - 2 * S
  return(Q_0 ^ exp(unprop))
}

# Define integrand 2
integrand2 <- function(X2, S) {
  unprop <- 0.5 + 0.7 * X2 - 2 * S
  return(Q_0 ^ exp(unprop))
}

# Compute integrals using pracma::integral2
result1 <- integral2(integrand1, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
result2 <- integral2(integrand2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)

# Final expectation
expectation <- 0.5 * result1$Q + 0.5 * result2$Q
print(expectation)
result1








# # 定义参数
# k <- -0.001052
# 
# # 定义被积函数
# integrand <- function(u) {
#   exp(k * exp(u))
# }
# 
# # 计算积分
# result <- integrate(integrand, lower = 0, upper = 1.2)
# print(result$value)