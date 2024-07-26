library(mice)
library(dplyr)
library(ggplot2)
library(devtools)
library(SimEngine)
library(survival)
library(parallel)
library(truncnorm)



# This file simulates a data from a vaccine clinical trial, with biomarker(partially observed).
# It can be ran both locally and in a cluster.

# set up multi-cores
run_on_cluster(
  # use SimEngine
  first = {
    sim = new_sim()
    
    sim %<>% set_levels(
      n = 50000,
      surv_time = list(
        "Exp" = list(surv_type = "Exponential", surv_params = 1.5e-2), # may be some problems
        "Gom" = list(surv_type = "Gompertz", surv_params = c(0.2138, 7e-8))
      ),
      error_meth = list(
        "Add" = list(meth_type = "additive", meth_params = c(0, 359.1)),
        "Mul" = list(meth_type = "multiplicative", meth_params = c(-0.2029, 0.637))
      ),
      error_relation = c("Classical", "Berkson")
    )
    
    create_data = function(n, meth_type, meth_params, error_relation, surv_type, surv_params) {
      # X1
      X1 = rnorm(n = n, mean = 24.3, sd = 8.38)
      
      # X2
      X2_true = rnorm(n = n, mean = 266.84, sd = 507.82)
      if (error_relation == "Classical") {
        # errors
        if (meth_type == "additive") {
          error = rnorm(n = n, mean = meth_params[1], sd = meth_params[2])
          X2 = X2_true + error
        } else if (meth_type == "multiplicative") {
          error = rlnorm(n = n, meanlog = meth_params[1], sdlog = meth_params[2])
          X2 = X2_true * error
        }
      } else if (error_relation == "Berkson") {
        # errors
        if (meth_type == "additive") {
          error = rnorm(n = n, mean = meth_params[1], sd = meth_params[2])
          X2 = X2_true - error
        } else if (meth_type == "multiplicative") {
          error = rlnorm(n = n, meanlog = meth_params[1], sdlog = meth_params[2])
          X2 = X2_true / error
        }
      }
      
      # biomarker
      S = rtruncnorm(n = n, a = 0, b = 1, mean = 0.5, sd = 0.2) # truncated normal in (0, 1)
      prob_tmp = 1 / (1 + exp(- 0.05*X1 + 0.005*X2)) # to make edge_prob distributed not too extremely
      val_tmp = rbinom(n, prob = prob_tmp, size = 1)
      S = (1 - val_tmp) * S
      
      # survival time
      U = runif(n = n)
      if (surv_type == "Exponential") {
        lambda = surv_params
        t = -log(U) / (lambda * exp(0.15*X1 + 0.001*X2_true + 5*S))
      } else if (surv_type == "Gompertz") {
        alpha = surv_params[1]
        lambda = surv_params[2]
        t = 1/alpha * log(1 - (alpha * log(U)) / (lambda * exp(0.15*X1 + 0.001*X2_true + 5*S)))
      }
      
      # censore time
      U = runif(n = n)
      if (surv_type == "Exponential") {
        lambda = surv_params
        C = - log(U) / (lambda * exp(0.15*X1 + 0.001*X2_true))
      } else if (surv_type == "Gompertz") {
        alpha = surv_params[1]
        lambda = surv_params[2]
        C = (- 1/lambda * log(U) * exp(0.15*X1 + 0.001*X2_true)) ^ (alpha)
      }
      
      # delta
      delta = ifelse(t <= C, 1, 0)
      
      # observed time
      Y = pmin(t, C)
      
      # two-phase indicator
      t0 = 200 # set the time of interest
      prob_tmp = 1 / (1 + exp(0.15*X1 + 0.001*X2_true - 1))
      Z = delta*I(Y <= t0) + (1 - delta*I(Y <= t0)) * rbinom(n = n, size = 1, prob = prob_tmp)
      
      data = data.frame("Y" = Y, "delta" = delta, "S" = S, "X1" = X1, "X2" = X2, "Z" = Z)
      return(data)
    }
    
    sim %<>% set_config(num_sim = 1000, n_cores = 4, seed = 1018)
    
    sim %<>% set_script(function() {
      library(survival)
      dat_phaseOne = create_data(L$n, L$error_meth$meth_type, L$error_meth$meth_params, L$error_relation, L$surv_time$surv_type, L$surv_time$surv_params)
      dat_phaseTwo = dat_phaseOne %>%
        filter(Z == 1) # use phase two data
      model = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_phaseTwo)
      return(list(
        "beta_X1_hat" = model$coef["X1"],
        "beta_X2_hat" = model$coef["X2"],
        "beta_S_hat" = model$coef["S"],
        "beta_X1_se" = summary(model)$coef["X1", "se(coef)"],
        "beta_X2_se" = summary(model)$coef["X2", "se(coef)"],
        "beta_S_se" = summary(model)$coef["S", "se(coef)"],
        "X1_pctg" = (model$coef["X1"] - 0.15) / 0.15 * 100,
        "X2_pctg" = (model$coef["X2"] - 0.001) / 0.001 *100,
        "S_pctg" = (model$coef["S"] - 5) / 5 *100,
        ".complex" = list(
          "model" = model,
          "cov_mtx" = vcov(model)
        )
      ))
    })
  },
  
  main = {
    sim %<>% run()
    print(sim$errors)
  },
  
  last = {
    # mean
    beta_X1 = sim %>% summarize(list(stat = "mean", x = "beta_X1_hat"))
    beta_X2 = sim %>% summarize(list(stat = "mean", x = "beta_X2_hat"))
    beta_S = sim %>% summarize(list(stat = "mean", x = "beta_S_hat"))
    # bias
    bias_beta = sim %>% summarize(
      list(stat = "bias", estimate = "beta_X1_hat", truth = 0.15, name = "bias_beta_X1"),
      list(stat = "bias", estimate = "beta_X2_hat", truth = 0.001, name = "bias_beta_X2"),
      list(stat = "bias", estimate = "beta_S_hat", truth = 0.001, name = "bias_beta_S"),
      list(stat = "coverage", estimate = "beta_X1_hat", se = "beta_X1_se", truth = 0.15, name = "cov_beta_X1"),
      list(stat = "coverage", estimate = "beta_X2_hat", se = "beta_X2_se", truth = 0.001, name = "cov_beta_X2"),
      list(stat = "coverage", estimate = "beta_S_hat", se = "beta_S_se", truth = 5, name = "cov_beta_S")
    )
    # bias percentage
    bias_beta_pct = sim %>% summarize(
      list(stat = "mean", x = "X1_pctg", name = "bias_beta_X1_pct"),
      list(stat = "mean", x = "X2_pctg", name = "bias_beta_X2_pct"),
      list(stat = "mean", x = "S_pctg", name = "bias_beta_S_pct")
    )
  },
  
  cluster_config = list(js = "slurm")
)
