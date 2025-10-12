# Two Phase Sampling estimated survival function
surv_two = function(model, t, data, type) {
  bh_func = basehaz(model, centered = F)
  index = which.min(abs(bh_func$time - t))
  H_0 = bh_func$hazard[index]
  beta = model$coefficients
  data$`I(S == 0)TRUE` = ifelse(data$S == 0, 1, 0)
  X_S = data[, c(names(model$coefficients))]
  unprop = exp(beta %*% t(X_S))
  Q = exp(- H_0 * unprop)
  
  if (type == "plc") {
    result = mean(Q)
  } else if (type == "vac") {
    result = sum(Q * data$ipw) / sum(data$ipw)
  } else if (type == "med") {
    X_S_med = X_S
    # X_S_med[, ncol(X_S_med)] = 0
    X_S_med$`S` = 0
    X_S_med$`I(S == 0)TRUE` = 1
    unprop_med = exp(beta %*% t(X_S_med))
    Q_med = exp(- H_0 * unprop_med)
    result = sum(Q_med * data$ipw) / sum(data$ipw)
  }
  return(result)
}