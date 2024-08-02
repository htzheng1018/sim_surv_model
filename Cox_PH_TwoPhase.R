library(SimEngine)



# This file simulates a data with a partially observed biomarker.
# Two-phase sampling method
# It can be ran both locally and in a cluster.



# start time
start_time = Sys.time()

# set up multi-cores
run_on_cluster(
  # use SimEngine
  first = {
    sim = new_sim()
    
    sim %<>% set_levels(
      n = c(100, 200, 500, 1000),
      surv_time = list(
        "Exp" = list(surv_type = "Exponential", surv_params = 1.5e-2), # may be some problems
        "Gom" = list(surv_type = "Gompertz", surv_params = c(0.2138, 7e-8))
      )
    )
    
    create_data = function(n, surv_type, surv_params) {
      # id
      id = seq(1, n)
      
      # treatment
      treat = sample(0:1, n, replace = TRUE, prob = c(0.7, 0.3))
      
      # X1
      X1 = rnorm(n = n, mean = 24.3, sd = 8.38)
      
      # X2
      X2 = rnorm(n = n, mean = 266.84, sd = 507.82)
      
      # biomarker
      S = rtruncnorm(n = n, a = 0, b = 1, mean = 0.5, sd = 0.2) # truncated normal in (0, 1)
      prob_tmp = 1 / (1 + exp(- 0.05*X1 + 0.005*X2 + treat)) # to make edge_prob distributed not too extremely
      val_tmp = rbinom(n, prob = prob_tmp, size = 1)
      S = (1 - val_tmp) * S
      
      # survival time
      U = runif(n = n)
      if (surv_type == "Exponential") {
        lambda = surv_params
        t = -log(U) / (lambda * exp(0.15*X1 + 0.001*X2 - 5*S + treat))
      } else if (surv_type == "Gompertz") {
        alpha = surv_params[1]
        lambda = surv_params[2]
        t = 1/alpha * log(1 - (alpha * log(U)) / (lambda * exp(0.15*X1 + 0.001*X2 - 5*S + treat)))
      }
      
      # censore time
      U = runif(n = n)
      if (surv_type == "Exponential") {
        lambda = surv_params
        C = - log(U) / (lambda * exp(0.15*X1 + 0.001*X2))
      } else if (surv_type == "Gompertz") {
        alpha = surv_params[1]
        lambda = surv_params[2]
        C = (- 1/lambda * log(U) * exp(0.15*X1 + 0.001*X2)) ^ (alpha)
      }
      
      # delta
      delta = ifelse(t <= C, 1, 0)
      
      # observed time
      Y = pmin(t, C)
      
      # two-phase indicator
      t0 = 200 # set the time of interest
      prob_tmp = 1 / (1 + exp(0.15*X1 + 0.001*X2 - 1))
      Z = delta*I(Y <= t0) + (1 - delta*I(Y <= t0)) * rbinom(n = n, size = 1, prob = prob_tmp)
      
      # temporary dataframe
      data = data.frame("id" = id, "treat" = treat, "Y" = Y, "delta" = delta, "S" = S, "X1" = X1, "X2" = X2, "Z" = Z)
      
      # using ipwpoint function to generate inverse probability weights
      ip_weights = ipwpoint(
        exposure = Z,
        family = "binomial",  # The treatment is binary
        link = "logit",
        denominator = ~ treat + delta + S + X1 + X2,
        data = data
      )$ipw.weights
      
      # final data
      data = data %>% 
        dplyr::mutate(ipw = ip_weights)
      
      return(data)
    }
    
    sim %<>% set_config(num_sim = 1000, n_cores = 4, seed = 1018,
                        packages = c("survival", "parallel", "truncnorm", "devtools", "ipw")
                        )
    
    sim %<>% set_script(function() {
      dat_phaseOne = create_data(L$n, L$surv_time$surv_type, L$surv_time$surv_params)
      dat_phaseTwo = dat_phaseOne %>%
        dplyr::filter(Z == 1) # use phase two data
      model = coxph(Surv(Y, delta) ~ X1 + X2 + S + treat, data = dat_phaseTwo, weights = ipw)
      return(list(
        "beta_X1_hat" = model$coef["X1"],
        "beta_X2_hat" = model$coef["X2"],
        "beta_S_hat" = model$coef["S"],
        "beta_treat_hat" = model$coef["treat"],
        "beta_X1_se" = summary(model)$coef["X1", "se(coef)"],
        "beta_X2_se" = summary(model)$coef["X2", "se(coef)"],
        "beta_S_se" = summary(model)$coef["S", "se(coef)"],
        "beta_treat_se" = summary(model)$coef["treat", "se(coef)"],
        "X1_pctg" = (model$coef["X1"] - 0.15) / 0.15 * 100,
        "X2_pctg" = (model$coef["X2"] - 0.001) / 0.001 *100,
        "S_pctg" = (model$coef["S"] + 5) / (-5) *100,
        "treat_pctg" = (model$coef["treat"] - 1) / 1 *100,
        ".complex" = list(
          "model" = model,
          "cov_mtx" = vcov(model),
          "data" = dat_phaseTwo
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
    beta_X1 = sim %>% SimEngine::summarize(list(stat = "mean", x = "beta_X1_hat"))
    beta_X2 = sim %>% SimEngine::summarize(list(stat = "mean", x = "beta_X2_hat"))
    beta_S = sim %>% SimEngine::summarize(list(stat = "mean", x = "beta_S_hat"))
    beta_treat = sim %>% SimEngine::summarize(list(stat = "mean", x = "beta_treat_hat"))
    # bias
    bias_beta = sim %>% SimEngine::summarize(
      list(stat = "bias", estimate = "beta_X1_hat", truth = 0.15, name = "bias_beta_X1"),
      list(stat = "bias", estimate = "beta_X2_hat", truth = 0.001, name = "bias_beta_X2"),
      list(stat = "bias", estimate = "beta_S_hat", truth = -5, name = "bias_beta_S"),
      list(stat = "bias", estimate = "beta_treat_hat", truth = 1, name = "bias_beta_treat"),
      list(stat = "coverage", estimate = "beta_X1_hat", se = "beta_X1_se", truth = 0.15, name = "cov_beta_X1"),
      list(stat = "coverage", estimate = "beta_X2_hat", se = "beta_X2_se", truth = 0.001, name = "cov_beta_X2"),
      list(stat = "coverage", estimate = "beta_S_hat", se = "beta_S_se", truth = -5, name = "cov_beta_S"),
      list(stat = "coverage", estimate = "beta_treat_hat", se = "beta_treat_se", truth = 1, name = "cov_beta_treat")
    )
    # bias percentage
    bias_beta_pct = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "X1_pctg", name = "bias_beta_X1_pct"),
      list(stat = "mean", x = "X2_pctg", name = "bias_beta_X2_pct"),
      list(stat = "mean", x = "S_pctg", name = "bias_beta_S_pct"),
      list(stat = "mean", x = "treat_pctg", name = "bias_beta_treat_pct")
    )
  },
  
  cluster_config = list(js = "slurm")
)



# end time
end_time = Sys.time()
execution_time = end_time - start_time
print(execution_time)
