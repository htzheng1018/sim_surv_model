



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
      n = c(500, 1000),
      # n = c(500, 1000, 2000, 4000, 8000),
      surv_time = list(
        "Exp" = list(surv_type = "Exponential", surv_params = 2e-2),
        "Gom" = list(surv_type = "Gompertz", surv_params = c(0.1, 1e-3))
      )
    )
    
    sim %<>% set_config(num_sim = 1000, n_cores = 4, seed = 1018,
                        packages = c("survival", "parallel", "truncnorm", "devtools", "ipw")
    )
    
    sim %<>% set_script(function() {
      dat_phaseOne = create_data(L$n, L$surv_time$surv_type, L$surv_time$surv_params, "iid") # phase one data (original)
      dat_phaseTwo_vac = dat_phaseOne %>%
        dplyr::filter(Z == 1 & treat==1) # vaccine group in phase two data
      dat_phaseOne_plc = dat_phaseOne[dat_phaseOne$treat == 0, ] # treat = 0 in placebo group
      dat_phaseOne_vac = dat_phaseOne[dat_phaseOne$treat == 1, ] # treat = 1 in vaccine group
      model_two_plc = coxph(Surv(Y, delta) ~ X1 + X2, data = dat_phaseOne_plc) # no s in placebo group
      model_two_vac = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_phaseTwo_vac, weights = ipw) # s in vaccine group
      
      
      
      # choose a specific time
      # time_max = round(max(dat_phaseOne$Y))
      # true_plc = c()
      # true_vac = c()
      # for (i in 1: time_max) {
      #   true_plc[i] = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, i, dat_phaseOne, "plc")
      #   true_vac[i] = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, i, dat_phaseOne_vac, "vac")
      # }
      # t_plc = which.min(abs(true_plc - 0.5))
      # print("placebo:")
      # print(t_plc)
      # t_vac = which.min(abs(true_vac - 0.5))
      # print("vaccine:")
      # print(t_vac)
      
      if (L$surv_time$surv_type == "Exponential") {
        t_plc = 19
        t_vac = 42
      } else if (L$surv_time$surv_type == "Gompertz") {
        t_plc = 36
        t_vac = 44
      }
      
      
      
      # bootstrap to get the variance of true survival functions and estimators
      surv_ci_plc = boot_ci(dat_phaseOne_plc, t_plc, "plc") # variance in placebo group
      surv_ci_vac = boot_ci(dat_phaseTwo_vac, t_vac, "vac") # variance in vaccine group
      
      # get the Survival probability at the specific time point
      Q_true_plc = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, t_plc, dat_phaseOne, "plc", "math")
      Q_true_vac = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, t_vac, dat_phaseOne_vac, "vac", "math")
      Q_est_km_plc = surv_km(t_plc, dat_phaseOne_plc) # km estimator for placebo group
      Q_est_km_vac = surv_km(t_vac, dat_phaseOne_vac) # km estimator for vaccine group
      Q_est_two_plc = surv_two(model_two_plc, t_plc, dat_phaseOne_plc)
      Q_est_two_vac = surv_two(model_two_vac, t_vac, dat_phaseTwo_vac)
      
      # get the true SE
      se_est_km = se_km(t, dat_phaseOne)
      se_est_two = se_two(t, dat_phaseOne)
      
      return(list(
        # survival functions
        "Q_true_plc" = Q_true_plc,
        "Q_true_vac" = Q_true_vac,
        "Q_est_km_plc" = Q_est_km_plc,
        "Q_est_km_vac" = Q_est_km_vac,
        "Q_est_two_plc" = Q_est_two_plc,
        "Q_est_two_vac" = Q_est_two_vac,
        "km_pctg_plc" = (Q_est_km_plc - Q_true_plc) / Q_true_plc * 100,
        "km_pctg_vac" = (Q_est_km_vac - Q_true_vac) / Q_true_vac * 100,
        "two_pctg_plc" = (Q_est_two_plc - Q_true_plc) / Q_true_plc * 100,
        "two_pctg_vac" = (Q_est_two_vac - Q_true_vac) / Q_true_vac * 100,
        
        # variance in placebo group
        "km_low_plc" = surv_ci_plc$km_low,
        "km_up_plc" = surv_ci_plc$km_up,
        "se_km_plc_boot" = surv_ci_plc$km_se,
        "se_km_vac_boot" = surv_ci_vac$km_se,
        # "se_km_est" = se_est_km,
        "two_low_plc" = surv_ci_plc$two_low,
        "two_up_plc" = surv_ci_plc$two_up,
        "se_two_plc_boot" = surv_ci_plc$two_se,
        "se_two_vac_boot" = surv_ci_vac$two_se,
        # "se_two_est" = se_est_two,
        
        # variance in vaccine group
        "km_low_vac" = surv_ci_vac$km_low,
        "km_up_vac" = surv_ci_vac$km_up,
        "se_km_vac_boot" = surv_ci_vac$km_se,
        # "se_km_est" = se_est_km,
        "two_low_vac" = surv_ci_vac$two_low,
        "two_up_vac" = surv_ci_vac$two_up,
        "se_two_vac_boot" = surv_ci_vac$two_se,
        # "se_two_est" = se_est_two,
        
        
        # "se_km_pctg" = (surv_ci$km_se - se_est_km) / se_est_km * 100,
        # "se_two_pctg" = (surv_ci$two_se - se_est_two) / se_est_two * 100,
        ".complex" = list(
          "model_plc" = model_two_plc,
          "data_plc" = dat_phaseOne_plc,
          "data_vac" = dat_phaseTwo_vac
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
    Q_true = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "Q_true_plc"),
      list(stat = "mean", x = "Q_true_vac")
    )
    Q_est_km = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "Q_est_km_plc"),
      list(stat = "mean", x = "Q_est_km_vac")
    )
    Q_est_two = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "Q_est_two_plc"),
      list(stat = "mean", x = "Q_est_two_vac")
    )
    
    # bias
    bias_Q_km = sim %>% SimEngine::summarize(
      list(stat = "bias", estimate = "Q_est_km_plc", truth = "Q_true_plc", name = "bias_km_plc"),
      list(stat = "bias", estimate = "Q_est_km_vac", truth = "Q_true_vac", name = "bias_km_vac")
    )
    bias_Q_two = sim %>% SimEngine::summarize(
      list(stat = "bias", estimate = "Q_est_two_plc", truth = "Q_true_plc", name = "bias_twophase_plc"),
      list(stat = "bias", estimate = "Q_est_two_vac", truth = "Q_true_vac", name = "bias_twophase_vac")
    )
    # bias percentage
    bias_Q_pct_km = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "km_pctg_plc", name = "bias_km_pct_plc"),
      list(stat = "mean", x = "km_pctg_vac", name = "bias_km_pct_vac")
    )
    bias_Q_pct_two = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "two_pctg_plc", name = "bias_twophase_pct_plc"),
      list(stat = "mean", x = "two_pctg_vac", name = "bias_twophase_pct_vac")
    )
    
    # SE accuracy
    accuracy_se_km = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "se_km_plc_boot", name = "se_km_plc"),
      list(stat = "mean", x = "se_km_vac_boot", name = "se_km_vac")
      # list(stat = "mean", x = "se_km_pctg", name = "se_bias_km_pct"),
      # list(stat = "mean", x = "se_km_pctg", name = "se_bias_km_pct")
    )
    accuracy_se_two = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "se_two_plc_boot", name = "se_two_plc"),
      list(stat = "mean", x = "se_two_vac_boot", name = "se_two_vac")
      # list(stat = "mean", x = "se_two_pctg", name = "se_bias_two_pct"),
      # list(stat = "mean", x = "se_two_pctg", name = "se_bias_two_pct")
    )
    
    # coverage
    coverage_km = sim %>% SimEngine::summarize(
      list(stat = "coverage", lower = "km_low_plc", upper = "km_up_plc", truth = "Q_true_plc", name = "cov_km_plc"),
      list(stat = "coverage", lower = "km_low_vac", upper = "km_up_vac", truth = "Q_true_vac", name = "cov_km_vac")
    )
    coverage_two = sim %>% SimEngine::summarize(
      list(stat = "coverage", lower = "two_low_plc", upper = "two_up_plc", truth = "Q_true_plc", name = "cov_twophase_plc"),
      list(stat = "coverage", lower = "two_low_vac", upper = "two_up_vac", truth = "Q_true_vac", name = "cov_twophase_vac")
    )
  },
  
  cluster_config = list(js = "slurm")
)



# end time
end_time = Sys.time()
execution_time = end_time - start_time
print(execution_time)




