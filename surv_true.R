# true survival function (placebo)
surv_true_plc = function(surv_type, surv_params, t, data) {
  if (surv_type == "Exponential") {
    lambda = surv_params
    Q_0 = exp(- lambda * t)
  } else if ( surv_type == "Gompertz") {
    alpha = surv_params[1]
    lambda = surv_params[2]
    Q_0 = exp(lambda/alpha * (1 - exp(alpha*t)))
  }
  # unprop = 0.15 * data$X1 + 0.001 * data$X2 - 5 * data$S + data$treat
  # unprop = 0.15 * data$X1 + 0.001 * data$X2 - 5 * data$S
  unprop = 0.15 * data$X1 + 0.001 * data$X2
  result = Q_0 ^ (exp(unprop))
  return(mean(result))
}

# true survival function (vaccine)
surv_true_vac = function(surv_type, surv_params, t, data) {
  if (surv_type == "Exponential") {
    lambda = surv_params
    Q_0 = exp(- lambda * t)
  } else if ( surv_type == "Gompertz") {
    alpha = surv_params[1]
    lambda = surv_params[2]
    Q_0 = exp(lambda/alpha * (1 - exp(alpha*t)))
  }
  # unprop = 0.15 * data$X1 + 0.001 * data$X2 - 5 * data$S + data$treat
  unprop = 0.15 * data$X1 + 0.001 * data$X2 - 5 * data$S
  # unprop = 0.15 * data$X1 + 0.001 * data$X2
  result = Q_0 ^ (exp(unprop))
  return(mean(result))
}




surv_0 = function(surv_type, surv_params, t, data) {
  if (surv_type == "Exponential") {
    lambda = surv_params
    Q_0 = exp(- lambda * t)
  } else if ( surv_type == "Gompertz") {
    alpha = surv_params[1]
    lambda = surv_params[2]
    Q_0 = exp(lambda/alpha * (1 - exp(alpha*t)))
  }
  return(Q_0)
}