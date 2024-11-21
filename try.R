


################
#### set up ####
################

# load SimEngine + functions
{
  library(SimEngine)
  source("create_data.R", local = T)
  source("surv_true.R", local = T)
  source("surv_km.R", local = T)
  source("surv_two.R", local = T)
  source("se_km.R", local = T)
  source("se_two.R", local = T)
  source("boot_ci.R", local = T)
}


##############
#### MAIN ####
##############

# It can be ran both locally and in a cluster.

# Fisrt step, for survival function in placebo group, verify the bias and the coverage of it using a Bootstrap method.

# start time
start_time = Sys.time()

# set up multi-cores
run_on_cluster(
  # use SimEngine
  first = {
    sim = new_sim()
    
    sim %<>% set_levels(
      # n = c(1000, 2000, 4000, 8000),
      n = c(500, 1000, 2000),
      surv_time = list(
        "Exp" = list(surv_type = "Exponential", surv_params = 1.5e-2), # may be some problems
        "Gom" = list(surv_type = "Gompertz", surv_params = c(0.2138, 7e-8))
      )
    )
    
    sim %<>% set_config(num_sim = 100, n_cores = 4, seed = 1018,
                        packages = c("survival", "parallel", "truncnorm", "devtools", "ipw")
    )
    
    sim %<>% set_script(function() {
      dat_phaseOne = create_data(L$n, L$surv_time$surv_type, L$surv_time$surv_params)
      model_one = coxph(Surv(Y, delta) ~ X1 + X2 + treat, data = dat_phaseOne)
      
      dat_phaseTwo = dat_phaseOne %>%
        dplyr::filter(Z == 1) # use phase two data
      dat_two_plc = dat_phaseTwo[dat_phaseTwo$treat == 0, ] # treat = 0 in placebo group
      model_two = coxph(Surv(Y, delta) ~ X1 + X2, data = dat_two_plc, weights = ipw) # no s in placebo group
      
      # choose a specific time
      time_max = round(max(dat_phaseOne$Y))
      true = c()
      for (i in 1: time_max) {
        true[i] = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, i, dat_phaseOne)
      }
      t = which.min(abs(true - 0.75))
      
      # bootstrap to get the variance of true survival functions and estimators
      surv_ci = boot_ci(dat_two_plc, t)
      
      # get the Survival probability at the specific time point
      Q_true = surv_0(L$surv_time$surv_type, L$surv_time$surv_params, t, dat_phaseOne)
      Q_est_two = surv_two(model_two, t, dat_two_plc)
      
      # get the true SE
      # se_est_km = se_km(t, dat_phaseOne)
      # se_est_two = se_two(t, dat_phaseOne)
      
      return(list(
        "Q_true" = Q_true,
        "Q_est_two" = Q_est_two,
        "two_low" = surv_ci$two_low,
        "two_up" = surv_ci$two_up,
        # "se_two_boot" = surv_ci$two_se,
        # "se_two_est" = se_est_two,
        # "two_pctg" = (Q_est_two - Q_true) / Q_true * 100,
        # "se_two_pctg" = (surv_ci$two_se - se_est_two) / se_est_two * 100,
        ".complex" = list(
          "model" = model_two,
          "data" = dat_two_plc
          # "ci" = surv_ci
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
    Q_true = sim %>% SimEngine::summarize(list(stat = "mean", x = "Q_true"))
    Q_est_two = sim %>% SimEngine::summarize(list(stat = "mean", x = "Q_est_two"))
    # bias
    bias_Q = sim %>% SimEngine::summarize(
      list(stat = "bias", estimate = "Q_est_two", truth = "Q_true", name = "bias_twophase"),
      list(stat = "coverage", lower = "two_low", upper = "two_up", truth = "Q_true", name = "cov_twophase")
    )
    # bias percentage
    # bias_Q_pct = sim %>% SimEngine::summarize(
    #   list(stat = "mean", x = "two_pctg", name = "bias_twophase_pct")
    # )
    # SE accuracy
    # accuracy_se = sim %>% SimEngine::summarize(
    #   list(stat = "mean", x = "se_two_boot", name = "se_two"),
    #   list(stat = "mean", x = "se_two_pctg", name = "se_bias_two_pct")
    # )
  },
  
  cluster_config = list(js = "slurm")
)



# end time
end_time = Sys.time()
execution_time = end_time - start_time
print(execution_time)







