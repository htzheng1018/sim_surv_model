library(mice)
library(dplyr)
library(ggplot2)
library(devtools)
library(SimEngine)
library(survival)
library(parallel)



# This file simulates data of Cox PH model using multi-cores.
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
      # age
      age = rnorm(n = n, mean = 24.3, sd = 8.38)
      
      # radon
      radon_true = rnorm(n = n, mean = 266.84, sd = 507.82)
      if (error_relation == "Classical") {
        # errors
        if (meth_type == "additive") {
          error = rnorm(n = n, mean = meth_params[1], sd = meth_params[2])
          radon = radon_true + error
        } else if (meth_type == "multiplicative") {
          error = rlnorm(n = n, meanlog = meth_params[1], sdlog = meth_params[2])
          radon = radon_true * error
        }
      } else if (error_relation == "Berkson") {
        # errors
        if (meth_type == "additive") {
          error = rnorm(n = n, mean = meth_params[1], sd = meth_params[2])
          radon = radon_true - error
        } else if (meth_type == "multiplicative") {
          error = rlnorm(n = n, meanlog = meth_params[1], sdlog = meth_params[2])
          radon = radon_true / error
        }
      }
    
      # uniform
      U = runif(n = n)
      
      # survival time
      if (surv_type == "Exponential") {
        lambda = surv_params
        t = -log(U) / (lambda * exp(0.15*age + 0.001*radon_true))
      } else if (surv_type == "Gompertz") {
        alpha = surv_params[1]
        lambda = surv_params[2]
        t = 1/alpha * log(1 - (alpha * log(U)) / (lambda * exp(0.15*age + 0.001*radon_true)))
      }
      
      # delta
      delta = rep(1, n)
      
      data = data.frame("X" = t, "delta" = delta, "radon" = radon, "age" = age)
      return(data)
    }

    sim %<>% set_config(num_sim = 1000, n_cores = 4, seed = 1018)
    
    sim %<>% set_script(function() {
      library(survival)
      dat = create_data(L$n, L$error_meth$meth_type, L$error_meth$meth_params, L$error_relation, L$surv_time$surv_type, L$surv_time$surv_params)
      model = coxph(Surv(X, delta) ~ radon + age, data = dat)
      return(list(
        "beta_radon_hat" = model$coef["radon"],
        "beta_age_hat" = model$coef["age"],
        "beta_radon_se" = summary(model)$coef["radon", "se(coef)"],
        "beta_age_se" = summary(model)$coef["age", "se(coef)"],
        "radon_pctg" = (model$coef["radon"] - 0.001) / 0.001 *100,
        "age_pctg" = (model$coef["age"] - 0.15) / 0.15 * 100,
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
    beta_radon = sim %>% summarize(list(stat = "mean", x = "beta_radon_hat"))
    beta_age = sim %>% summarize(list(stat = "mean", x = "beta_age_hat"))
    # bias
    bias_beta = sim %>% summarize(
      list(stat = "bias", estimate = "beta_radon_hat", truth = 0.001, name = "bias_beta_radon"),
      list(stat = "bias", estimate = "beta_age_hat", truth = 0.15, name = "bias_beta_age"),
      list(stat = "coverage", estimate = "beta_radon_hat", se = "beta_radon_se", truth = 0.001, name = "cov_beta_radon"),
      list(stat = "coverage", estimate = "beta_age_hat", se = "beta_age_se", truth = 0.15, name = "cov_beta_age")
    )
    # bias percentage
    bias_beta_pct = sim %>% summarize(
      list(stat = "mean", x = "radon_pctg", name = "bias_beta_radon_pct"),
      list(stat = "mean", x = "age_pctg", name = "bias_beta_age_pct")
    )
  },
  
  cluster_config = list(js = "slurm")
)



                     