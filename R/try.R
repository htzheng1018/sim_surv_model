# This could be executed both locally and on a computing cluster.

# For the Two-phase sampling estimators, we assessed their bias and variance within a CoxPH framework.  
# Additionally, we used a Bootstrap method to evaluate the coverage and variance bias of the estimators.



################
#### set up ####
################

# load SimEngine + functions
{
  library(SimEngine)
  library(kableExtra)
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
      n = c(500, 1000),
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
      # normal data and normal model, without the indicator
      dat_phaseOne = create_data(L$n, L$surv_time$surv_type, L$surv_time$surv_params, "complex") # phase one data (original)
      dat_phaseTwo_vac = dat_phaseOne %>%
        dplyr::filter(Z == 1 & treat==1) # vaccine group in phase two data
      dat_phaseOne_plc = dat_phaseOne[dat_phaseOne$treat == 0, ] # treat = 0 in placebo group
      dat_phaseOne_vac = dat_phaseOne[dat_phaseOne$treat == 1, ] # treat = 1 in vaccine group
      model_two_plc = coxph(Surv(Y, delta) ~ X1 + X2, data = dat_phaseOne_plc) # no s in placebo group
      model_two_vac = coxph(Surv(Y, delta) ~ X1 + X2 + S, data = dat_phaseTwo_vac, weights = ipw) # s in vaccine group
      
      # normal data with biomarker indicator model
      model_two_vac_plus = coxph(Surv(Y, delta) ~ X1 + X2 + S + I(S == 0), data = dat_phaseTwo_vac, weights = ipw) # refined model
      
      # biomarker indicator data and biomarker indicator model
      dat_phaseOne_pro = create_data(L$n, L$surv_time$surv_type, L$surv_time$surv_params, "complex", ind = T) # phase one data
      dat_phaseTwo_vac_pro = dat_phaseOne_pro %>%
        dplyr::filter(Z == 1 & treat==1) # vaccine group in phase two data
      dat_phaseOne_plc_pro = dat_phaseOne_pro[dat_phaseOne_pro$treat == 0, ] # treat = 0 in placebo group
      dat_phaseOne_vac_pro = dat_phaseOne_pro[dat_phaseOne_pro$treat == 1, ] # treat = 1 in vaccine group
      model_two_plc_pro = coxph(Surv(Y, delta) ~ X1 + X2, data = dat_phaseOne_plc_pro) # no s in placebo group
      model_two_vac_pro = coxph(Surv(Y, delta) ~ X1 + X2 + S + I(S == 0), data = dat_phaseTwo_vac_pro, weights = ipw) # refined model
      
      
      
      # # choose a specific time
      # time_max = round(max(dat_phaseOne$Y))
      # true_plc = c()
      # true_vac = c()
      # true_med = c()
      # for (i in 1: time_max) {
      #   true_plc[i] = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, i, dat_phaseOne_ind, "plc", "sample", ind = T)
      #   true_vac[i] = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, i, dat_phaseOne_ind_vac, "vac", "sample", ind = T)
      #   true_med[i] = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, i, dat_phaseOne_ind_vac, "med", "sample", ind = T)
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
        t_plc_pro = 19
        t_vac_pro = 46
        t_med_pro = 32
      } else if (L$surv_time$surv_type == "Gompertz") {
        t_plc = 37
        t_vac = 44
        t_med = 37
        t_plc_pro = 37
        t_vac_pro = 45
        t_med_pro = 42
      }
      
      
      
      # bootstrap to get the variance of true survival functions and estimators
      # variance in normal data with normal model
      ci_boot_plc = ci(dat_phaseOne_plc, t_plc, "plc", "bootstrap")
      ci_boot_vac = ci(dat_phaseTwo_vac, t_vac, "vac", "bootstrap")
      ci_boot_med = ci(dat_phaseTwo_vac, t_med, "med", "bootstrap")
      
      # variance in normal data with biomarker indicator model
      ci_boot_plc_plus = ci(dat_phaseOne_plc, t_plc, "plc", "bootstrap", ind = T)
      ci_boot_vac_plus = ci(dat_phaseTwo_vac, t_vac, "vac", "bootstrap", ind = T)
      ci_boot_med_plus = ci(dat_phaseTwo_vac, t_med, "med", "bootstrap", ind = T)
      
      # variance in biomarker indicator data with biomarker indicator model
      ci_boot_plc_pro = ci(dat_phaseOne_plc_pro, t_plc_pro, "plc", "bootstrap", ind = T)
      ci_boot_vac_pro = ci(dat_phaseTwo_vac_pro, t_vac_pro, "vac", "bootstrap", ind = T)
      ci_boot_med_pro = ci(dat_phaseTwo_vac_pro, t_med_pro, "med", "bootstrap", ind = T)
      
      
      
      # get the Survival probability at the specific time point
      # Kaplan-Meier functions
      Q_est_km_plc = surv_km(t_plc, dat_phaseOne_plc, "plc")
      Q_est_km_vac = surv_km(t_vac, dat_phaseTwo_vac, "vac")
      Q_est_km_med = surv_km(t_med, dat_phaseTwo_vac, "med")
      Q_est_km_plc_pro = surv_km(t_plc, dat_phaseOne_plc_pro, "plc")
      Q_est_km_vac_pro = surv_km(t_vac, dat_phaseTwo_vac_pro, "vac")
      Q_est_km_med_pro = surv_km(t_med, dat_phaseTwo_vac_pro, "med")
      
      # true survival functions
      Q_true_plc = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, t_plc, dat_phaseOne, "plc", "math")
      Q_true_vac = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, t_vac, dat_phaseOne_vac, "vac", "math")
      Q_true_med = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, t_med, dat_phaseOne_vac, "med", "math")
      Q_true_plc_pro = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, t_plc_pro, dat_phaseOne_pro, "plc", "math", ind = T)
      Q_true_vac_pro = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, t_vac_pro, dat_phaseOne_vac_pro, "vac", "math", ind = T)
      Q_true_med_pro = surv_true(L$surv_time$surv_type, L$surv_time$surv_params, t_med_pro, dat_phaseOne_vac_pro, "med", "math", ind = T)
      
      # survival functions in normal data with normal model
      Q_est_two_plc = surv_two(model_two_plc, t_plc, dat_phaseOne_plc, "plc")
      Q_est_two_vac = surv_two(model_two_vac, t_vac, dat_phaseTwo_vac, "vac")
      Q_est_two_med = surv_two(model_two_vac, t_med, dat_phaseTwo_vac, "med")
      
      # survival functions in normal data with biomarker indicator model
      Q_est_two_plc_plus = surv_two(model_two_plc, t_plc, dat_phaseOne_plc, "plc", ind = T)
      Q_est_two_vac_plus = surv_two(model_two_vac_plus, t_vac, dat_phaseTwo_vac, "vac", ind = T)
      Q_est_two_med_plus = surv_two(model_two_vac_plus, t_med, dat_phaseTwo_vac, "med", ind = T)
      
      # survival functions in biomarker indicator data with biomarker indicator model
      Q_est_two_plc_pro = surv_two(model_two_plc_pro, t_plc_pro, dat_phaseOne_plc_pro, "plc", ind = T)
      Q_est_two_vac_pro = surv_two(model_two_vac_pro, t_vac_pro, dat_phaseTwo_vac_pro, "vac", ind = T)
      Q_est_two_med_pro = surv_two(model_two_vac_pro, t_med_pro, dat_phaseTwo_vac_pro, "med", ind = T)
      
      # get the true SE
      # se_est_km = se_km(t, dat_phaseOne)
      # se_est_two = se_two(t, dat_phaseOne)
      
      
      
      return(list(
        # true survival functions
        "Q_true_plc" = Q_true_plc,
        "Q_true_vac" = Q_true_vac,
        "Q_true_med" = Q_true_med,
        "Q_true_plc_pro" = Q_true_plc_pro,
        "Q_true_vac_pro" = Q_true_vac_pro,
        "Q_true_med_pro" = Q_true_med_pro,
        
        # estimators
        # Kaplan-Meier estimator
        "Q_est_km_plc" = Q_est_km_plc,
        "Q_est_km_vac" = Q_est_km_vac,
        "Q_est_km_med" = Q_est_km_med,
        "Q_est_km_plc_pro" = Q_est_km_plc_pro,
        "Q_est_km_vac_pro" = Q_est_km_vac_pro,
        "Q_est_km_med_pro" = Q_est_km_med_pro,
        # two-phase sampling estimators
        "Q_est_two_plc" = Q_est_two_plc,
        "Q_est_two_vac" = Q_est_two_vac,
        "Q_est_two_med" = Q_est_two_med,
        "Q_est_two_plc_plus" = Q_est_two_plc_plus,
        "Q_est_two_vac_plus" = Q_est_two_vac_plus,
        "Q_est_two_med_plus" = Q_est_two_med_plus,
        "Q_est_two_plc_pro" = Q_est_two_plc_pro,
        "Q_est_two_vac_pro" = Q_est_two_vac_pro,
        "Q_est_two_med_pro" = Q_est_two_med_pro,
        
        # bias
        # K-M
        "km_pctg_plc" = (Q_est_km_plc - Q_true_plc) / Q_true_plc * 100,
        "km_pctg_vac" = (Q_est_km_vac - Q_true_vac) / Q_true_vac * 100,
        "km_pctg_med" = (Q_est_km_med - Q_true_med) / Q_true_med * 100,
        "km_pctg_plc_pro" = (Q_est_km_plc_pro - Q_true_plc_pro) / Q_true_plc_pro * 100,
        "km_pctg_vac_pro" = (Q_est_km_vac_pro - Q_true_vac_pro) / Q_true_vac_pro * 100,
        "km_pctg_med_pro" = (Q_est_km_med_pro - Q_true_med_pro) / Q_true_med_pro * 100,
        # two-phase sampling
        "two_pctg_plc" = (Q_est_two_plc - Q_true_plc) / Q_true_plc * 100,
        "two_pctg_vac" = (Q_est_two_vac - Q_true_vac) / Q_true_vac * 100,
        "two_pctg_med" = (Q_est_two_med - Q_true_med) / Q_true_med * 100,
        "two_pctg_plc_plus" = (Q_est_two_plc_plus - Q_true_plc) / Q_true_plc * 100,
        "two_pctg_vac_plus" = (Q_est_two_vac_plus - Q_true_vac) / Q_true_vac * 100,
        "two_pctg_med_plus" = (Q_est_two_med_plus - Q_true_med) / Q_true_med * 100,
        "two_pctg_plc_pro" = (Q_est_two_plc_pro - Q_true_plc_pro) / Q_true_plc_pro * 100,
        "two_pctg_vac_pro" = (Q_est_two_vac_pro - Q_true_vac_pro) / Q_true_vac_pro * 100,
        "two_pctg_med_pro" = (Q_est_two_med_pro - Q_true_med_pro) / Q_true_med_pro * 100,
        
        # variance
        # K-M
        "km_plc_low" = ci_boot_plc$km_low,
        "km_plc_up" = ci_boot_plc$km_up,
        "se_km_plc_boot" = ci_boot_plc$km_se,
        "var_km_plc_boot" = ci_boot_plc$km_se ^ 2,
        "km_vac_low" = ci_boot_vac$km_low,
        "km_vac_up" = ci_boot_vac$km_up,
        "se_km_vac_boot" = ci_boot_vac$km_se,
        "var_km_vac_boot" = ci_boot_vac$km_se ^ 2,
        "km_plc_low_pro" = ci_boot_plc_pro$km_low,
        "km_plc_up_pro" = ci_boot_plc_pro$km_up,
        "se_km_plc_boot_pro" = ci_boot_plc_pro$km_se,
        "var_km_plc_boot_pro" = ci_boot_plc_pro$km_se ^ 2,
        "km_vac_low_pro" = ci_boot_vac_pro$km_low,
        "km_vac_up_pro" = ci_boot_vac_pro$km_up,
        "se_km_vac_boot_pro" = ci_boot_vac_pro$km_se,
        "var_km_vac_boot_pro" = ci_boot_vac_pro$km_se ^ 2,
        # two-phase sampling
        "two_plc_low" = ci_boot_plc$two_low,
        "two_plc_up" = ci_boot_plc$two_up,
        "se_two_plc_boot" = ci_boot_plc$two_se,
        "var_two_plc_boot" = ci_boot_plc$two_se ^ 2,
        "two_plc_low_plus" = ci_boot_plc_plus$two_low,
        "two_plc_up_plus" = ci_boot_plc_plus$two_up,
        "se_two_plc_boot_plus" = ci_boot_plc_plus$two_se,
        "var_two_plc_boot_plus" = ci_boot_plc_plus$two_se ^ 2,
        "two_plc_low_pro" = ci_boot_plc_pro$two_low,
        "two_plc_up_pro" = ci_boot_plc_pro$two_up,
        "se_two_plc_boot_pro" = ci_boot_plc_pro$two_se,
        "var_two_plc_boot_pro" = ci_boot_plc_pro$two_se ^ 2,
        "two_vac_low" = ci_boot_vac$two_low,
        "two_vac_up" = ci_boot_vac$two_up,
        "se_two_vac_boot" = ci_boot_vac$two_se,
        "var_two_vac_boot" = ci_boot_vac$two_se ^ 2,
        "two_vac_low_plus" = ci_boot_vac_plus$two_low,
        "two_vac_up_plus" = ci_boot_vac_plus$two_up,
        "se_two_vac_boot_plus" = ci_boot_vac_plus$two_se,
        "var_two_vac_boot_plus" = ci_boot_vac_plus$two_se ^ 2,
        "two_vac_low_pro" = ci_boot_vac_pro$two_low,
        "two_vac_up_pro" = ci_boot_vac_pro$two_up,
        "se_two_vac_boot_pro" = ci_boot_vac_pro$two_se,
        "var_two_vac_boot_pro" = ci_boot_vac_pro$two_se ^ 2,
        "two_med_low" = ci_boot_med$two_low,
        "two_med_up" = ci_boot_med$two_up,
        "se_two_med_boot" = ci_boot_med$two_se,
        "var_two_med_boot" = ci_boot_med$two_se ^ 2,
        "two_med_low_plus" = ci_boot_med_plus$two_low,
        "two_med_up_plus" = ci_boot_med_plus$two_up,
        "se_two_med_boot_plus" = ci_boot_med_plus$two_se,
        "var_two_med_boot_plus" = ci_boot_med_plus$two_se ^ 2,
        "two_med_low_pro" = ci_boot_med_pro$two_low,
        "two_med_up_pro" = ci_boot_med_pro$two_up,
        "se_two_med_boot_pro" = ci_boot_med_pro$two_se,
        "var_two_med_boot_pro" = ci_boot_med_pro$two_se ^ 2,
        
        # "se_km_est" = se_est_km,
        # "se_two_est" = se_est_two,
        # "se_km_est" = se_est_km,
        # "se_two_est" = se_est_two,
        # "se_km_est" = se_est_km,
        # "se_two_est" = se_est_two,
        
        ".complex" = list(
          "model_plc" = model_two_plc,
          "model_plc_pro" = model_two_plc_pro,
          "model_vac" = model_two_vac,
          "model_vac_plus" = model_two_vac_plus,
          "model_vac_pro" = model_two_vac_pro,
          "data_plc" = dat_phaseOne_plc,
          "data_plc_pro" = dat_phaseOne_plc_pro,
          "data_vac" = dat_phaseTwo_vac,
          "data_vac_pro" = dat_phaseTwo_vac_pro
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
      list(stat = "mean", x = "Q_true_med"),
      list(stat = "mean", x = "Q_true_plc_pro"),
      list(stat = "mean", x = "Q_true_vac_pro"),
      list(stat = "mean", x = "Q_true_med_pro")
    )
    Q_est_km = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "Q_est_km_plc"),
      list(stat = "mean", x = "Q_est_km_vac"),
      list(stat = "mean", x = "Q_est_km_med"),
      list(stat = "mean", x = "Q_est_km_plc_pro"),
      list(stat = "mean", x = "Q_est_km_vac_pro"),
      list(stat = "mean", x = "Q_est_km_med_pro")
    )
    Q_est_two = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "Q_est_two_plc"),
      list(stat = "mean", x = "Q_est_two_vac"),
      list(stat = "mean", x = "Q_est_two_med"),
      list(stat = "mean", x = "Q_est_two_plc_plus"),
      list(stat = "mean", x = "Q_est_two_vac_plus"),
      list(stat = "mean", x = "Q_est_two_med_plus"),
      list(stat = "mean", x = "Q_est_two_plc_pro"),
      list(stat = "mean", x = "Q_est_two_vac_pro"),
      list(stat = "mean", x = "Q_est_two_med_pro")
    )
    
    # bias
    bias_Q_km = sim %>% SimEngine::summarize(
      list(stat = "bias", estimate = "Q_est_km_plc", truth = "Q_true_plc", name = "bias_km_plc"),
      list(stat = "bias", estimate = "Q_est_km_vac", truth = "Q_true_vac", name = "bias_km_vac"),
      list(stat = "bias", estimate = "Q_est_km_med", truth = "Q_true_med", name = "bias_km_med"),
      list(stat = "bias", estimate = "Q_est_km_plc_pro", truth = "Q_true_plc_pro", name = "bias_km_plc_pro"),
      list(stat = "bias", estimate = "Q_est_km_vac_pro", truth = "Q_true_vac_pro", name = "bias_km_vac_pro"),
      list(stat = "bias", estimate = "Q_est_km_med_pro", truth = "Q_true_med_pro", name = "bias_km_med_pro")
    )
    bias_Q_two = sim %>% SimEngine::summarize(
      list(stat = "bias", estimate = "Q_est_two_plc", truth = "Q_true_plc", name = "bias_twophase_plc"),
      list(stat = "bias", estimate = "Q_est_two_vac", truth = "Q_true_vac", name = "bias_twophase_vac"),
      list(stat = "bias", estimate = "Q_est_two_med", truth = "Q_true_med", name = "bias_twophase_med"),
      list(stat = "bias", estimate = "Q_est_two_plc_plus", truth = "Q_true_plc", name = "bias_twophase_plc_plus"),
      list(stat = "bias", estimate = "Q_est_two_vac_plus", truth = "Q_true_vac", name = "bias_twophase_vac_plus"),
      list(stat = "bias", estimate = "Q_est_two_med_plus", truth = "Q_true_med", name = "bias_twophase_med_plus"),
      list(stat = "bias", estimate = "Q_est_two_plc_pro", truth = "Q_true_plc_pro", name = "bias_twophase_plc_pro"),
      list(stat = "bias", estimate = "Q_est_two_vac_pro", truth = "Q_true_vac_pro", name = "bias_twophase_vac_pro"),
      list(stat = "bias", estimate = "Q_est_two_med_pro", truth = "Q_true_med_pro", name = "bias_twophase_med_pro")
    )
    # bias percentage
    bias_Q_pct_km = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "km_pctg_plc", name = "bias_km_pct_plc"),
      list(stat = "mean", x = "km_pctg_vac", name = "bias_km_pct_vac"),
      list(stat = "mean", x = "km_pctg_med", name = "bias_km_pct_med"),
      list(stat = "mean", x = "km_pctg_plc_pro", name = "bias_km_pct_plc_pro"),
      list(stat = "mean", x = "km_pctg_vac_pro", name = "bias_km_pct_vac_pro"),
      list(stat = "mean", x = "km_pctg_med_pro", name = "bias_km_pct_med_pro")
    )
    bias_Q_pct_two = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "two_pctg_plc", name = "bias_twophase_pct_plc"),
      list(stat = "mean", x = "two_pctg_vac", name = "bias_twophase_pct_vac"),
      list(stat = "mean", x = "two_pctg_med", name = "bias_twophase_pct_med"),
      list(stat = "mean", x = "two_pctg_plc_plus", name = "bias_twophase_pct_plc_plus"),
      list(stat = "mean", x = "two_pctg_vac_plus", name = "bias_twophase_pct_vac_plus"),
      list(stat = "mean", x = "two_pctg_med_plus", name = "bias_twophase_pct_med_plus"),
      list(stat = "mean", x = "two_pctg_plc_pro", name = "bias_twophase_pct_plc_pro"),
      list(stat = "mean", x = "two_pctg_vac_pro", name = "bias_twophase_pct_vac_pro"),
      list(stat = "mean", x = "two_pctg_med_pro", name = "bias_twophase_pct_med_pro")
    )
    
    # variance in average
    var_km = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "var_km_plc_boot", name = "var_km_plc_boot"),
      list(stat = "mean", x = "var_km_vac_boot", name = "var_km_vac_boot"),
      list(stat = "mean", x = "var_km_plc_boot_pro", name = "var_km_plc_boot_pro"),
      list(stat = "mean", x = "var_km_vac_boot_pro", name = "var_km_vac_boot_pro")
    )
    var_two = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "var_two_plc_boot", name = "var_two_plc_boot"),
      list(stat = "mean", x = "var_two_vac_boot", name = "var_two_vac_boot"),
      list(stat = "mean", x = "var_two_med_boot", name = "var_two_med_boot"),
      list(stat = "mean", x = "var_two_plc_boot_plus", name = "var_two_plc_boot_plus"),
      list(stat = "mean", x = "var_two_vac_boot_plus", name = "var_two_vac_boot_plus"),
      list(stat = "mean", x = "var_two_med_boot_plus", name = "var_two_med_boot_plus"),
      list(stat = "mean", x = "var_two_plc_boot_pro", name = "var_two_plc_boot_pro"),
      list(stat = "mean", x = "var_two_vac_boot_pro", name = "var_two_vac_boot_pro"),
      list(stat = "mean", x = "var_two_med_boot_pro", name = "var_two_med_boot_pro")
    )
    
    # coverage
    coverage_km = sim %>% SimEngine::summarize(
      list(stat = "coverage", lower = "km_plc_low", upper = "km_plc_up", truth = "Q_true_plc", name = "cov_km_plc"),
      list(stat = "coverage", lower = "km_vac_low", upper = "km_vac_up", truth = "Q_true_vac", name = "cov_km_vac"),
      list(stat = "coverage", lower = "km_plc_low_pro", upper = "km_plc_up_pro", truth = "Q_true_plc_pro", name = "cov_km_plc_pro"),
      list(stat = "coverage", lower = "km_vac_low_pro", upper = "km_vac_up_pro", truth = "Q_true_vac_pro", name = "cov_km_vac_pro")
    )
    coverage_two = sim %>% SimEngine::summarize(
      list(stat = "coverage", lower = "two_plc_low", upper = "two_plc_up", truth = "Q_true_plc", name = "cov_twophase_plc"),
      list(stat = "coverage", lower = "two_vac_low", upper = "two_vac_up", truth = "Q_true_vac", name = "cov_twophase_vac"),
      list(stat = "coverage", lower = "two_med_low", upper = "two_med_up", truth = "Q_true_med", name = "cov_twophase_med"),
      list(stat = "coverage", lower = "two_plc_low_plus", upper = "two_plc_up_plus", truth = "Q_true_plc", name = "cov_twophase_plc_plus"),
      list(stat = "coverage", lower = "two_vac_low_plus", upper = "two_vac_up_plus", truth = "Q_true_vac", name = "cov_twophase_vac_plus"),
      list(stat = "coverage", lower = "two_med_low_plus", upper = "two_med_up_plus", truth = "Q_true_med", name = "cov_twophase_med_plus"),
      list(stat = "coverage", lower = "two_plc_low_pro", upper = "two_plc_up_pro", truth = "Q_true_plc_pro", name = "cov_twophase_plc_pro"),
      list(stat = "coverage", lower = "two_vac_low_pro", upper = "two_vac_up_pro", truth = "Q_true_vac_pro", name = "cov_twophase_vac_pro"),
      list(stat = "coverage", lower = "two_med_low_pro", upper = "two_med_up_pro", truth = "Q_true_med_pro", name = "cov_twophase_med_pro")
    )
  },
  
  cluster_config = list(js = "slurm")
)



# end time
end_time = Sys.time()
execution_time = end_time - start_time
print(execution_time)





