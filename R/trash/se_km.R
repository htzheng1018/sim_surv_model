# SE for Kaplan Meier
se_km = function(t, data) {
  km.est = survfit(Surv(Y, delta) ~ 1, data = data, conf.type = "log-log")
  index = which.min(abs(km.est$time - t))
  result = km.est$std.err[index]
  return(result)
}
