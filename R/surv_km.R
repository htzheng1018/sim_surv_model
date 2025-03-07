# Kaplan Meier estimated survival function
surv_km = function(t, data, type) {
  if (type == "plc") {
    km.est = survfit(Surv(Y, delta) ~ 1, data = data, conf.type = "log-log")
  } else if (type == "vac") {
    km.est = survfit(Surv(Y, delta) ~ 1, data = data, conf.type = "log-log", weights = data$ipw)
  }
  index = which.min(abs(km.est$time - t))
  result = km.est$surv[index]
  return(result)
}