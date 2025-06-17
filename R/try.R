# This could be executed both locally and on a computing cluster.

# For the Two-phase sampling estimators, we assessed their bias and variance within a CoxPH framework.  
# Additionally, we used a Bootstrap method to evaluate the coverage and variance bias of the estimators.



################
#### set up ####
################

# load SimEngine + functions
{
  library(SimEngine)
  source("R/create_data.R", local = T)
  source("R/surv_true.R", local = T)
  source("R/surv_km.R", local = T)
  source("R/surv_two.R", local = T)
  source("R/se_km.R", local = T)
  source("R/se_two.R", local = T)
  source("R/ci.R", local = T)
}



##############
#### MAIN ####
##############

# start time
start_time = Sys.time()

# set up multi-cores
run_on_cluster(
  # use SimEngine
  first = {
    sim = new_sim()
    
    sim %<>% set_levels(
      n = 1000,
      # n = c(500, 1000, 2000, 4000, 8000),
      surv_time = list(
        "Exp" = list(surv_type = "Exponential", surv_params = 2e-2),
        "Gom" = list(surv_type = "Gompertz", surv_params = c(0.1, 1e-3))
      )
    )
    
    sim %<>% set_config(num_sim = 1000, n_cores = 4, seed = 1018,
                        packages = c("survival", "parallel", "truncnorm", "devtools", "ipw", "pracma")
    )
    
    sim %<>% set_script(function() {
      dat_phaseOne = create_data(L$n, L$surv_time$surv_type, L$surv_time$surv_params, "complex") # phase one data (original)
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
      # true_med = c()
      # for (i in 1: time_max) {
      #   true_plc[i] = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, i, dat_phaseOne, "plc", "sample")
      #   true_vac[i] = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, i, dat_phaseOne_vac, "vac", "sample")
      #   true_med[i] = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, i, dat_phaseOne_vac, "med", "sample")
      # }
      # t_plc = which.min(abs(true_plc - 0.5))
      # print("placebo:")
      # print(t_plc)
      # t_vac = which.min(abs(true_vac - 0.5))
      # print("vaccine:")
      # print(t_vac)
      # t_med = which.min(abs(true_med - 0.5))
      # print("mediation:")
      # print(t_med)
      
      if (L$surv_time$surv_type == "Exponential") {
        t_plc = 19
        t_vac = 42
        t_med = 19
      } else if (L$surv_time$surv_type == "Gompertz") {
        t_plc = 37
        t_vac = 44
        t_med = 37
      }
      
      
      
      # bootstrap to get the variance of true survival functions and estimators
      ci_boot_plc = ci(dat_phaseOne_plc, t_plc, "plc", "bootstrap") # variance in placebo group
      ci_boot_vac = ci(dat_phaseTwo_vac, t_vac, "vac", "bootstrap") # variance in vaccine group
      ci_boot_med = ci(dat_phaseTwo_vac, t_med, "med", "bootstrap") # variance in vaccine group
      
      # get the Survival probability at the specific time point
      Q_true_plc = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, t_plc, dat_phaseOne, "plc", "math")
      Q_true_vac = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, t_vac, dat_phaseOne_vac, "vac", "math")
      Q_true_med = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, t_med, dat_phaseOne_vac, "med", "math")
      Q_est_km_plc = surv_km(t_plc, dat_phaseOne_plc, "plc") # km estimator for placebo group
      Q_est_km_vac = surv_km(t_vac, dat_phaseTwo_vac, "vac") # km estimator for vaccine group
      Q_est_km_med = surv_km(t_med, dat_phaseTwo_vac, "med") # km estimator for mediation group
      Q_est_two_plc = surv_two(model_two_plc, t_plc, dat_phaseOne_plc, "plc")
      Q_est_two_vac = surv_two(model_two_vac, t_vac, dat_phaseTwo_vac, "vac")
      Q_est_two_med = surv_two(model_two_vac, t_med, dat_phaseTwo_vac, "med")
      
      # get the true SE
      # se_est_km = se_km(t, dat_phaseOne)
      # se_est_two = se_two(t, dat_phaseOne)
      
      
      
      return(list(
        # survival functions
        "Q_true_plc" = Q_true_plc,
        "Q_true_vac" = Q_true_vac,
        "Q_true_med" = Q_true_med,
        "Q_est_km_plc" = Q_est_km_plc,
        "Q_est_km_vac" = Q_est_km_vac,
        "Q_est_km_med" = Q_est_km_med,
        "Q_est_two_plc" = Q_est_two_plc,
        "Q_est_two_vac" = Q_est_two_vac,
        "Q_est_two_med" = Q_est_two_med,
        "km_pctg_plc" = (Q_est_km_plc - Q_true_plc) / Q_true_plc * 100,
        "km_pctg_vac" = (Q_est_km_vac - Q_true_vac) / Q_true_vac * 100,
        "km_pctg_med" = (Q_est_km_med - Q_true_med) / Q_true_med * 100,
        "two_pctg_plc" = (Q_est_two_plc - Q_true_plc) / Q_true_plc * 100,
        "two_pctg_vac" = (Q_est_two_vac - Q_true_vac) / Q_true_vac * 100,
        "two_pctg_med" = (Q_est_two_med - Q_true_med) / Q_true_med * 100,
        
        # variance in placebo group
        "km_plc_low" = ci_boot_plc$km_low,
        "km_plc_up" = ci_boot_plc$km_up,
        "se_km_plc_boot" = ci_boot_plc$km_se,
        # "se_km_est" = se_est_km,
        "two_plc_low" = ci_boot_plc$two_low,
        "two_plc_up" = ci_boot_plc$two_up,
        "se_two_plc_boot" = ci_boot_plc$two_se,
        # "se_two_est" = se_est_two,
        
        # variance in vaccine group
        "km_vac_low" = ci_boot_vac$km_low,
        "km_vac_up" = ci_boot_vac$km_up,
        "se_km_vac_boot" = ci_boot_vac$km_se,
        # "se_km_est" = se_est_km,
        "two_vac_low" = ci_boot_vac$two_low,
        "two_vac_up" = ci_boot_vac$two_up,
        "se_two_vac_boot" = ci_boot_vac$two_se,
        # "se_two_est" = se_est_two,
        
        # variance in mediation group
        "km_med_low" = ci_boot_med$km_low,
        "km_med_up" = ci_boot_med$km_up,
        "se_km_med_boot" = ci_boot_med$km_se,
        # "se_km_est" = se_est_km,
        "two_med_low" = ci_boot_med$two_low,
        "two_med_up" = ci_boot_med$two_up,
        "se_two_med_boot" = ci_boot_med$two_se,
        # "se_two_est" = se_est_two,
        
        ".complex" = list(
          "model_plc" = model_two_plc,
          "model_vac" = model_two_vac,
          "model_med" = model_two_vac,
          "data_plc" = dat_phaseOne_plc,
          "data_vac" = dat_phaseTwo_vac,
          "data_med" = dat_phaseTwo_vac
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
      list(stat = "mean", x = "Q_true_vac"),
      list(stat = "mean", x = "Q_true_med")
    )
    Q_est_km = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "Q_est_km_plc"),
      list(stat = "mean", x = "Q_est_km_vac"),
      list(stat = "mean", x = "Q_est_km_med")
    )
    Q_est_two = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "Q_est_two_plc"),
      list(stat = "mean", x = "Q_est_two_vac"),
      list(stat = "mean", x = "Q_est_two_med")
    )
    
    # bias
    bias_Q_km = sim %>% SimEngine::summarize(
      list(stat = "bias", estimate = "Q_est_km_plc", truth = "Q_true_plc", name = "bias_km_plc"),
      list(stat = "bias", estimate = "Q_est_km_vac", truth = "Q_true_vac", name = "bias_km_vac"),
      list(stat = "bias", estimate = "Q_est_km_med", truth = "Q_true_med", name = "bias_km_med")
    )
    bias_Q_two = sim %>% SimEngine::summarize(
      list(stat = "bias", estimate = "Q_est_two_plc", truth = "Q_true_plc", name = "bias_twophase_plc"),
      list(stat = "bias", estimate = "Q_est_two_vac", truth = "Q_true_vac", name = "bias_twophase_vac"),
      list(stat = "bias", estimate = "Q_est_two_med", truth = "Q_true_med", name = "bias_twophase_med")
    )
    # bias percentage
    bias_Q_pct_km = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "km_pctg_plc", name = "bias_km_pct_plc"),
      list(stat = "mean", x = "km_pctg_vac", name = "bias_km_pct_vac"),
      list(stat = "mean", x = "km_pctg_med", name = "bias_km_pct_med")
    )
    bias_Q_pct_two = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "two_pctg_plc", name = "bias_twophase_pct_plc"),
      list(stat = "mean", x = "two_pctg_vac", name = "bias_twophase_pct_vac"),
      list(stat = "mean", x = "two_pctg_med", name = "bias_twophase_pct_med")
    )
    
    # SE accuracy
    # accuracy_se_km = sim %>% SimEngine::summarize(
    #   list(stat = "mean", x = "se_km_plc_boot", name = "se_km_plc"),
    #   list(stat = "mean", x = "se_km_vac_boot", name = "se_km_vac"),
    # list(stat = "mean", x = "se_km_pctg", name = "se_bias_km_pct"),
    # list(stat = "mean", x = "se_km_pctg", name = "se_bias_km_pct")
    # )
    # accuracy_se_two = sim %>% SimEngine::summarize(
    #   list(stat = "mean", x = "se_two_plc_boot", name = "se_two_plc"),
    #   list(stat = "mean", x = "se_two_vac_boot", name = "se_two_vac"),
    # list(stat = "mean", x = "se_two_pctg", name = "se_bias_two_pct"),
    # list(stat = "mean", x = "se_two_pctg", name = "se_bias_two_pct")
    # )
    
    # coverage
    coverage_km = sim %>% SimEngine::summarize(
      list(stat = "coverage", lower = "km_plc_low", upper = "km_plc_up", truth = "Q_true_plc", name = "cov_km_plc"),
      list(stat = "coverage", lower = "km_vac_low", upper = "km_vac_up", truth = "Q_true_vac", name = "cov_km_vac"),
      list(stat = "coverage", lower = "km_med_low", upper = "km_med_up", truth = "Q_true_med", name = "cov_km_med")
    )
    coverage_two = sim %>% SimEngine::summarize(
      list(stat = "coverage", lower = "two_plc_low", upper = "two_plc_up", truth = "Q_true_plc", name = "cov_twophase_plc"),
      list(stat = "coverage", lower = "two_vac_low", upper = "two_vac_up", truth = "Q_true_vac", name = "cov_twophase_vac"),
      list(stat = "coverage", lower = "two_med_low", upper = "two_med_up", truth = "Q_true_med", name = "cov_twophase_med")
    )
  },
  
  cluster_config = list(js = "slurm")
)



# save results
saveRDS(Q_est_km, file = "Evaluation/Q_est_km.rds")
saveRDS(Q_est_two, file = "Evaluation/Q_est_two.rds")
saveRDS(Q_true, file = "Evaluation/Q_true.rds")
saveRDS(bias_Q_km, file = "Evaluation/bias_Q_km.rds")
saveRDS(bias_Q_two, file = "Evaluation/bias_Q_two.rds")
saveRDS(bias_Q_pct_km, file = "Evaluation/bias_Q_pct_km.rds")
saveRDS(bias_Q_pct_two, file = "Evaluation/bias_Q_pct_two.rds")
saveRDS(coverage_km, file = "Evaluation/coverage_km.rds")
saveRDS(coverage_two, file = "Evaluation/coverage_two.rds")



# end time
end_time = Sys.time()
execution_time = end_time - start_time
print(execution_time)




