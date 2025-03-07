# SE for two-sphase sampling
se_two = function(t, data) {
  nn = nrow(data)
  R = 100
  surv_two.boot = c()
  
  for (r in 1: R) {
    boot.samp = sample(1: nn, size = nn, replace = TRUE)
    data.boot = data[boot.samp, ]
    model.boot = coxph(Surv(Y, delta) ~ X1 + X2 + treat, data = data.boot, weights = ipw)
    surv_two.boot[r] = surv_two(model.boot, t, data)
  }
  result = sqrt(sum((surv_two.boot - mean(surv_two.boot)) ^ 2) / (R - 1))
  return(result)
}