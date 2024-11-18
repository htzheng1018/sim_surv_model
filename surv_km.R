# Kaplan Meier estimated survival function
surv_km = function(t, data) {
  km.est = survfit(Surv(Y, delta) ~ 1, data = data, conf.type = "log-log")
  index = which.min(abs(km.est$time - t))
  result = km.est$surv[index]
  return(result)
}