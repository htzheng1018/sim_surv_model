set.seed(1018)
n = 1000
lambda = 2e-2
result11 = c()

for(i in (1: 180)) {
  Q_0 = exp(- lambda * i)
  integrand2 = function(X2, S) {
    unprop = 0.5*1 + 0.7*X2 - 2*S
    # Q = pmax(Q_0 ^ exp(unprop), 1e-2)
    Q = Q_0 ^ exp(unprop)
    prob_tmp = 1 / (1 + exp(0.5*1 + 0.7*X2 + 1))
    result = Q * (1 - prob_tmp) * dtruncnorm(S, a = 0, b = 1, mean = 0.5, sd = 0.2)
    return(result)
  }
  
  # result11[i] = integral2(function(X2, S) integrand2(X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
  # print(result11[i])
  
  tryCatch({
    result11[i] = integral2(function(X2, S) integrand2(X2, S), xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
  }, error = function(e) {result11[i] = result11[i - 1]}
  )
}








