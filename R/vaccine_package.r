
################
#### set up ####
################

# load SimEngine + functions
{
  library(SimEngine)
  library(kableExtra)
  source("R/create_data.R", local = T)
  source("R/true_func.r", local = T)
  source("R/est_tps.r", local = T)
  source("R/est_flx.r", local = T)
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
      # n = c(500, 1000, 2000, 4000, 8000),
      n = c(500, 1000),
      surv_time = list(
        "Exp" = list(surv_type = "Exponential", surv_params = 2e-2),
        "Gom" = list(surv_type = "Gompertz", surv_params = c(0.1, 1e-3))
      )
    )
    
    sim %<>% set_config(num_sim = 100, n_cores = 13, seed = 1018,
                        packages = c("survival", "parallel", "truncnorm", "devtools", "ipw", "pracma")
    )
    
    sim %<>% set_script(function() {
      # normal data and normal model, without the indicator
      dat_org = create_data(L$n, L$surv_time$surv_type, L$surv_time$surv_params, "complex") # phase one data (original)
      dat_ind = create_data(L$n, L$surv_time$surv_type, L$surv_time$surv_params, "complex", ind = T) # phase one data with biomarker indicator
      
      # choose a specific time
      if (L$surv_time$surv_type == "Exponential") {
        t_org = 30
        t_ind = 33
      } else if (L$surv_time$surv_type == "Gompertz") {
        t_org = 40
        t_ind = 42
      }
      
      # true value
      val_0_org = true_func(L$surv_time$surv_type, L$surv_time$surv_params, t_org, dat_org, "math")
      val_0_ind = true_func(L$surv_time$surv_type, L$surv_time$surv_params, t_ind, dat_ind, "math", ind = T)
      
      # two-phase sampling estimator and flexible estimator
      val_n_tps_org = est_tps(dat_org, t_org, boots = 100)
      val_n_tps_ind = est_tps(dat_ind, t_ind, boots = 100)
      val_n_flx_org = est_flx(dat_org, t_org, boots = 100)
      val_n_flx_ind = est_flx(dat_ind, t_ind, boots = 100)
      
      # results
      return(list(
        # true values
        "NIE_0_org" = val_0_org["NIE", "true"],
        "NIE_0_ind" = val_0_ind["NIE", "true"],
        "NDE_0_org" = val_0_org["NDE", "true"],
        "NDE_0_ind" = val_0_ind["NDE", "true"],
        "PM_0_org" = val_0_org["PM", "true"],
        "PM_0_ind" = val_0_ind["PM", "true"],

        # estimators
        "NIE_n_tps_org" = val_n_tps_org["NIE", "estimate"],
        "NIE_n_tps_ind" = val_n_tps_ind["NIE", "estimate"],
        "NDE_n_tps_org" = val_n_tps_org["NDE", "estimate"],
        "NDE_n_tps_ind" = val_n_tps_ind["NDE", "estimate"],
        "PM_n_tps_org" = val_n_tps_org["PM", "estimate"],
        "PM_n_tps_ind" = val_n_tps_ind["PM", "estimate"],
        "NIE_n_flx_org" = val_n_flx_org["NIE", "estimate"],
        "NIE_n_flx_ind" = val_n_flx_ind["NIE", "estimate"],
        "NDE_n_flx_org" = val_n_flx_org["NDE", "estimate"],
        "NDE_n_flx_ind" = val_n_flx_ind["NDE", "estimate"],
        "PM_n_flx_org" = val_n_flx_org["PM", "estimate"],
        "PM_n_flx_ind" = val_n_flx_ind["PM", "estimate"],
        
        # standard errors
        "NIE_se_tps_org" = val_n_tps_org["NIE", "se"],
        "NIE_se_tps_ind" = val_n_tps_ind["NIE", "se"],
        "NDE_se_tps_org" = val_n_tps_org["NDE", "se"],
        "NDE_se_tps_ind" = val_n_tps_ind["NDE", "se"],
        "PM_se_tps_org" = val_n_tps_org["PM", "se"],
        "PM_se_tps_ind" = val_n_tps_ind["PM", "se"],
        "NIE_se_flx_org" = val_n_flx_org["NIE", "se"],
        "NIE_se_flx_ind" = val_n_flx_ind["NIE", "se"],
        "NDE_se_flx_org" = val_n_flx_org["NDE", "se"],
        "NDE_se_flx_ind" = val_n_flx_ind["NDE", "se"],
        "PM_se_flx_org" = val_n_flx_org["PM", "se"],
        "PM_se_flx_ind" = val_n_flx_ind["PM", "se"],
        
        # bias percentage
        "NIE_bias_tps_org" = (val_n_tps_org["NIE", "estimate"] - val_0_org["NIE", "true"]) / val_0_org["NIE", "true"] * 100,
        "NIE_bias_tps_ind" = (val_n_tps_ind["NIE", "estimate"] - val_0_ind["NIE", "true"]) / val_0_ind["NIE", "true"] * 100,
        "NDE_bias_tps_org" = (val_n_tps_org["NDE", "estimate"] - val_0_org["NDE", "true"]) / val_0_org["NDE", "true"] * 100,
        "NDE_bias_tps_ind" = (val_n_tps_ind["NDE", "estimate"] - val_0_ind["NDE", "true"]) / val_0_ind["NDE", "true"] * 100,
        "PM_bias_tps_org" = (val_n_tps_org["PM", "estimate"] - val_0_org["PM", "true"]) / val_0_org["PM", "true"] * 100,
        "PM_bias_tps_ind" = (val_n_tps_ind["PM", "estimate"] - val_0_ind["PM", "true"]) / val_0_ind["PM", "true"] * 100,
        "NIE_bias_flx_org" = (val_n_flx_org["NIE", "se"] - val_0_org["NIE", "true"]) / val_0_org["NIE", "true"] * 100,
        "NIE_bias_flx_ind" = (val_n_flx_ind["NIE", "se"] - val_0_ind["NIE", "true"]) / val_0_ind["NIE", "true"] * 100,
        "NDE_bias_flx_org" = (val_n_flx_org["NDE", "se"] - val_0_org["NDE", "true"]) / val_0_org["NDE", "true"] * 100,
        "NDE_bias_flx_ind" = (val_n_flx_ind["NDE", "se"] - val_0_ind["NDE", "true"]) / val_0_ind["NDE", "true"] * 100,
        "PM_bias_flx_org" = (val_n_flx_org["PM", "se"] - val_0_org["PM", "true"]) / val_0_org["PM", "true"] * 100,
        "PM_bias_flx_ind" = (val_n_flx_ind["PM", "se"] - val_0_ind["PM", "true"]) / val_0_ind["PM", "true"] * 100,
        
        # 95% CI
        "NIE_low_tps_org" = val_n_tps_org["NIE", "low"],
        "NIE_low_tps_ind" = val_n_tps_ind["NIE", "low"],
        "NIE_up_tps_org" = val_n_tps_org["NIE", "up"],
        "NIE_up_tps_ind" = val_n_tps_ind["NIE", "up"],
        "NDE_low_tps_org" = val_n_tps_org["NDE", "low"],
        "NDE_low_tps_ind" = val_n_tps_ind["NDE", "low"],
        "NDE_up_tps_org" = val_n_tps_org["NDE", "up"],
        "NDE_up_tps_ind" = val_n_tps_ind["NDE", "up"],
        "PM_low_tps_org" = val_n_tps_org["PM", "low"],
        "PM_low_tps_ind" = val_n_tps_ind["PM", "low"],
        "PM_up_tps_org" = val_n_tps_org["PM", "up"],
        "PM_up_tps_ind" = val_n_tps_ind["PM", "up"],
        "NIE_low_flx_org" = val_n_flx_org["NIE", "low"],
        "NIE_low_flx_ind" = val_n_flx_ind["NIE", "low"],
        "NIE_up_flx_org" = val_n_flx_org["NIE", "up"],
        "NIE_up_flx_ind" = val_n_flx_ind["NIE", "up"],
        "NDE_low_flx_org" = val_n_flx_org["NDE", "low"],
        "NDE_low_flx_ind" = val_n_flx_ind["NDE", "low"],
        "NDE_up_flx_org" = val_n_flx_org["NDE", "up"],
        "NDE_up_flx_ind" = val_n_flx_ind["NDE", "up"],
        "PM_low_flx_org" = val_n_flx_org["PM", "low"],
        "PM_low_flx_ind" = val_n_flx_ind["PM", "low"],
        "PM_up_flx_org" = val_n_flx_org["PM", "up"],
        "PM_up_flx_ind" = val_n_flx_ind["PM", "up"],
        
        # other complex results
        ".complex" = list(
          "dat_org" = dat_org,
          "dat_ind" = dat_ind,
          "val_0_org" = val_0_org,
          "val_0_ind" = val_0_ind,
          "val_n_tps_org" = val_n_tps_org,
          "val_n_tps_ind" = val_n_tps_ind,
          "val_n_flx_org" = val_n_flx_org,
          "val_n_flx_ind" = val_n_flx_ind
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
    true_values = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "NIE_0_org"),
      list(stat = "mean", x = "NIE_0_ind"),
      list(stat = "mean", x = "NDE_0_org"),
      list(stat = "mean", x = "NDE_0_ind"),
      list(stat = "mean", x = "PM_0_org"),
      list(stat = "mean", x = "PM_0_ind")
    )
    estimators = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "NIE_n_tps_org"),
      list(stat = "mean", x = "NIE_n_tps_ind"),
      list(stat = "mean", x = "NDE_n_tps_org"),
      list(stat = "mean", x = "NDE_n_tps_ind"),
      list(stat = "mean", x = "PM_n_tps_org"),
      list(stat = "mean", x = "PM_n_tps_ind"),
      list(stat = "mean", x = "NIE_n_flx_org"),
      list(stat = "mean", x = "NIE_n_flx_ind"),
      list(stat = "mean", x = "NDE_n_flx_org"),
      list(stat = "mean", x = "NDE_n_flx_ind"),
      list(stat = "mean", x = "PM_n_flx_org"),
      list(stat = "mean", x = "PM_n_flx_ind")
    )
    
    # standard error
    standard_error = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "NIE_se_tps_org"),
      list(stat = "mean", x = "NIE_se_tps_ind"),
      list(stat = "mean", x = "NDE_se_tps_org"),
      list(stat = "mean", x = "NDE_se_tps_ind"),
      list(stat = "mean", x = "PM_se_tps_org"),
      list(stat = "mean", x = "PM_se_tps_ind"),
      list(stat = "mean", x = "NIE_se_flx_org"),
      list(stat = "mean", x = "NIE_se_flx_ind"),
      list(stat = "mean", x = "NDE_se_flx_org"),
      list(stat = "mean", x = "NDE_se_flx_ind"),
      list(stat = "mean", x = "PM_se_flx_org"),
      list(stat = "mean", x = "PM_se_flx_ind")
    )
    
    # bias
    bias = sim %>% SimEngine::summarize(
      list(stat = "bias", estimate = "NIE_n_tps_org", truth = "NIE_0_tps_org", name = "bias_NIE_tps_org"),
      list(stat = "bias", estimate = "NIE_n_tps_ind", truth = "NIE_0_tps_ind", name = "bias_NIE_tps_ind"),
      list(stat = "bias", estimate = "NDE_n_tps_org", truth = "NDE_0_tps_org", name = "bias_NDE_tps_org"),
      list(stat = "bias", estimate = "NDE_n_tps_ind", truth = "NDE_0_tps_ind", name = "bias_NDE_tps_ind"),
      list(stat = "bias", estimate = "PM_n_tps_org", truth = "PM_0_tps_org", name = "bias_PM_tps_org"),
      list(stat = "bias", estimate = "PM_n_tps_ind", truth = "PM_0_tps_ind", name = "bias_PM_tps_ind"),
      list(stat = "bias", estimate = "NIE_n_flx_org", truth = "NIE_0_flx_org", name = "bias_NIE_flx_org"),
      list(stat = "bias", estimate = "NIE_n_flx_ind", truth = "NIE_0_flx_ind", name = "bias_NIE_flx_ind"),
      list(stat = "bias", estimate = "NDE_n_flx_org", truth = "NDE_0_flx_org", name = "bias_NDE_flx_org"),
      list(stat = "bias", estimate = "NDE_n_flx_ind", truth = "NDE_0_flx_ind", name = "bias_NDE_flx_ind"),
      list(stat = "bias", estimate = "PM_n_flx_org", truth = "PM_0_flx_org", name = "bias_PM_flx_org"),
      list(stat = "bias", estimate = "PM_n_flx_ind", truth = "PM_0_flx_ind", name = "bias_PM_flx_ind")
    )
    
    # bias percentage
    bias_percentage = sim %>% SimEngine::summarize(
      list(stat = "mean", x = "NIE_bias_tps_org", name = "bias_pctg_NIE_tps_org"),
      list(stat = "mean", x = "NIE_bias_tps_ind", name = "bias_pctg_NIE_tps_ind"),
      list(stat = "mean", x = "NDE_bias_tps_org", name = "bias_pctg_NDE_tps_org"),
      list(stat = "mean", x = "NDE_bias_tps_ind", name = "bias_pctg_NDE_tps_ind"),
      list(stat = "mean", x = "PM_bias_tps_org", name = "bias_pctg_PM_tps_org"),
      list(stat = "mean", x = "PM_bias_tps_ind", name = "bias_pctg_PM_tps_ind"),
      list(stat = "mean", x = "NIE_bias_flx_org", name = "bias_pctg_NIE_flx_org"),
      list(stat = "mean", x = "NIE_bias_flx_ind", name = "bias_pctg_NIE_flx_ind"),
      list(stat = "mean", x = "NDE_bias_flx_org", name = "bias_pctg_NDE_flx_org"),
      list(stat = "mean", x = "NDE_bias_flx_ind", name = "bias_pctg_NDE_flx_ind"),
      list(stat = "mean", x = "PM_bias_flx_org", name = "bias_pctg_PM_flx_org"),
      list(stat = "mean", x = "PM_bias_flx_ind", name = "bias_pctg_PM_flx_ind")
    )
    
    # coverage
    coverage = sim %>% SimEngine::summarize(
      list(stat = "coverage", lower = "NIE_low_tps_org", upper = "NIE_up_tps_org", truth = "NIE_0_org", name = "cov_NIE_tps_org"),
      list(stat = "coverage", lower = "NIE_low_tps_ind", upper = "NIE_up_tps_ind", truth = "NIE_0_ind", name = "cov_NIE_tps_ind"),
      list(stat = "coverage", lower = "NDE_low_tps_org", upper = "NDE_up_tps_org", truth = "NDE_0_org", name = "cov_NDE_tps_org"),
      list(stat = "coverage", lower = "NDE_low_tps_ind", upper = "NDE_up_tps_ind", truth = "NDE_0_ind", name = "cov_NDE_tps_ind"),
      list(stat = "coverage", lower = "PM_low_tps_org", upper = "PM_up_tps_org", truth = "PM_0_org", name = "cov_PM_tps_org"),
      list(stat = "coverage", lower = "PM_low_tps_ind", upper = "PM_up_tps_ind", truth = "PM_0_ind", name = "cov_PM_tps_ind"),
      list(stat = "coverage", lower = "NIE_low_flx_org", upper = "NIE_up_flx_org", truth = "NIE_0_org", name = "cov_NIE_flx_org"),
      list(stat = "coverage", lower = "NIE_low_flx_ind", upper = "NIE_up_flx_ind", truth = "NIE_0_ind", name = "cov_NIE_flx_ind"),
      list(stat = "coverage", lower = "NDE_low_flx_org", upper = "NDE_up_flx_org", truth = "NDE_0_org", name = "cov_NDE_flx_org"),
      list(stat = "coverage", lower = "NDE_low_flx_ind", upper = "NDE_up_flx_ind", truth = "NDE_0_ind", name = "cov_NDE_flx_ind"),
      list(stat = "coverage", lower = "PM_low_flx_org", upper = "PM_up_flx_org", truth = "PM_0_org", name = "cov_PM_flx_org"),
      list(stat = "coverage", lower = "PM_low_flx_ind", upper = "PM_up_flx_ind", truth = "PM_0_ind", name = "cov_PM_flx_ind")
    )
  },
  
  cluster_config = list(js = "slurm")
)



# end time
end_time = Sys.time()
execution_time = end_time - start_time
print(execution_time)





