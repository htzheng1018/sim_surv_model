


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
      n = 1000,
      # n = c(500, 1000, 2000, 4000, 8000),
      surv_time = list(
        "Exp" = list(surv_type = "Exponential", surv_params = 1.5e-3), # may be some problems
        "Gom" = list(surv_type = "Gompertz", surv_params = c(0.2138, 7e-8))
      )
    )
    
    sim %<>% set_config(num_sim = 1000, n_cores = 4, seed = 1018,
                        packages = c("survival", "parallel", "truncnorm", "devtools", "ipw")
    )
    
    sim %<>% set_script(function() {
      dat_phaseOne = create_data(L$n, L$surv_time$surv_type, L$surv_time$surv_params, "complex")
      # model_one = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_phaseOne)
      
      dat_phaseTwo_vac = dat_phaseOne %>%
        dplyr::filter(Z == 1 & treat==1) # use phase two data
      dat_phaseOne_plc = dat_phaseOne[dat_phaseOne$treat == 0, ] # treat = 0 in placebo group
      dat_phaseOne_vac = dat_phaseOne[dat_phaseOne$treat == 1, ] # treat = 1 in vaccine group
      model_two_plc = coxph(Surv(Y, delta) ~ X1 + X2, data = dat_phaseOne_plc) # no s in placebo group
      model_two_vac = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_phaseTwo_vac, weights = ipw) # s in vaccine group
      
      # choose a specific time
      # time_max = round(max(dat_phaseOne$Y))
      # true = c()
      # for (i in 1: time_max) {
      #   true[i] = surv_true_plc(L$surv_time$surv_type, L$surv_time$surv_params, i, dat_phaseOne)
      # }
      # t = which.min(abs(true - 0.25))
      # print(t)
      if (L$surv_time$surv_type == "Exponential") {
        t = 10
      } else if (L$surv_time$surv_type == "Gompertz") {
        t = 50
      }
      
      # # bootstrap to get the variance of true survival functions and estimators
      # surv_ci = boot_ci(dat_two_plc, t)
      
      # get the Survival probability at the specific time point
      Q_true_plc = surv_true_plc(L$surv_time$surv_type, L$surv_time$surv_params, t, dat_phaseOne)
      # Q_true_vac = surv_true_vac(L$surv_time$surv_type, L$surv_time$surv_params, t, dat_phaseTwo_vac)
      Q_true_vac = surv_true_vac(L$surv_time$surv_type, L$surv_time$surv_params, t, dat_phaseOne_vac)
      Q_est_km = surv_km(t, dat_phaseOne_plc)
      Q_est_two_plc = surv_two(model_two_plc, t, dat_phaseOne_plc, "plc")
      Q_est_two_vac = surv_two(model_two_vac, t, dat_phaseTwo_vac, "vac")
      
      # get the true SE
      # se_est_km = se_km(t, dat_phaseOne)
      # se_est_two = se_two(t, dat_phaseOne)
      
      return(list(
        "Q_true_plc" = Q_true_plc,
        "Q_true_vac" = Q_true_vac,
        "Q_est_km" = Q_est_km,
        "Q_est_two_plc" = Q_est_two_plc,
        "Q_est_two_vac" = Q_est_two_vac,
        # "km_low" = surv_ci$km_low,
        # "km_up" = surv_ci$km_up,
        # "se_km_boot" = surv_ci$km_se,
        # "se_km_est" = se_est_km,
        # "two_low" = surv_ci$two_low,
        # "two_up" = surv_ci$two_up,
        # "se_two_boot" = surv_ci$two_se,
        # "se_two_est" = se_est_two,
        "km_pctg" = (Q_est_km - Q_true_plc) / Q_true_plc * 100,
        "two_pctg_plc" = (Q_est_two_plc - Q_true_plc) / Q_true_plc * 100,
        "two_pctg_vac" = (Q_est_two_vac - Q_true_vac) / Q_true_vac * 100,
        # "se_km_pctg" = (surv_ci$km_se - se_est_km) / se_est_km * 100,
        # "se_two_pctg" = (surv_ci$two_se - se_est_two) / se_est_two * 100,
        ".complex" = list(
          "model" = model_two_plc,
          "data" = dat_phaseOne_plc
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
    Q_true_plc = sim %>% SimEngine::summarize(list(stat = "mean", x = "Q_true_plc"))
    Q_true_vac = sim %>% SimEngine::summarize(list(stat = "mean", x = "Q_true_vac"))
    # Q_est_km = sim %>% SimEngine::summarize(list(stat = "mean", x = "Q_est_km"))
    Q_est_two_plc = sim %>% SimEngine::summarize(list(stat = "mean", x = "Q_est_two_plc"))
    Q_est_two_vac = sim %>% SimEngine::summarize(list(stat = "mean", x = "Q_est_two_vac"))
    # bias
    bias_Q = sim %>% SimEngine::summarize(
      # list(stat = "bias", estimate = "Q_est_km", truth = "Q_true", name = "bias_km"),
      list(stat = "bias", estimate = "Q_est_two_plc", truth = "Q_true_plc", name = "bias_twophase_plc"),
      list(stat = "bias", estimate = "Q_est_two_vac", truth = "Q_true_vac", name = "bias_twophase_vac")
      # list(stat = "coverage", lower = "km_low", upper = "km_up", truth = "Q_true", name = "cov_km"),
      # list(stat = "coverage", lower = "two_low", upper = "two_up", truth = "Q_true", name = "cov_twophase")
    )
    # bias percentage
    bias_Q_pct = sim %>% SimEngine::summarize(
      # list(stat = "mean", x = "km_pctg", name = "bias_km_pct"),
      list(stat = "mean", x = "two_pctg_plc", name = "bias_twophase_pct_plc"),
      list(stat = "mean", x = "two_pctg_vac", name = "bias_twophase_pct_vac")
    )
    # SE accuracy
    # accuracy_se = sim %>% SimEngine::summarize(
    #   list(stat = "mean", x = "se_km_boot", name = "se_km"),
    #   list(stat = "mean", x = "se_two_boot", name = "se_two"),
    #   list(stat = "mean", x = "se_km_pctg", name = "se_bias_km_pct"),
    #   list(stat = "mean", x = "se_two_pctg", name = "se_bias_two_pct")
    # )
    
    print(head(sim$results))
    print(bias_Q)
    print(bias_Q_pct)
    
  },
  
  cluster_config = list(js = "slurm")
)



# end time
end_time = Sys.time()
execution_time = end_time - start_time
print(execution_time)