
################
#### set up ####
################

# load SimEngine + functions
{
  library(SimEngine)
  library(vaccine)
  library(kableExtra)
  source("R/create_data.R", local = T)
  source("R/true_func.r", local = T)
  source("R/est_med.r", local = T)
}

# split the large result into the X-only and X+S versions
split_med_result = function(result, version = c("raw", "xs")) {
  version = match.arg(version)
  columns = if (version == "raw") {
    c("NIE_one", "NDE", "PM_one")
  } else if (version == "xs") {
    c("NIE_two", "NDE", "PM_two")
  }
  out = result[, columns, drop = FALSE]
  colnames(out) = c("NIE", "NDE", "PM")
  # if is not calculated for the X+S version
  if (version == "xs") {
    out[c("se_if", "low_if", "up_if"), ] = NA_real_
  }
  return(out)
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
      n = 8000,
      surv_time = list(
        "Exp" = list(surv_type = "Exponential", surv_params = 2e-2),
        "Gom" = list(surv_type = "Gompertz", surv_params = c(0.1, 1e-3))
      )
    )
    
    sim %<>% set_config(num_sim = 5, parallel = TRUE, n_cores = 13, seed = 1018,
                        packages = c("survival", "parallel", "truncnorm", "pracma", "dplyr", "vaccine")
    )
    
    sim %<>% set_script(function() {
      num_boot = 1000
      
      # normal data and normal model, without the indicator
      dat_org = create_data(L$n, L$surv_time$surv_type, L$surv_time$surv_params, "complex") # phase one data (original)
      dat_ind = create_data(L$n, L$surv_time$surv_type, L$surv_time$surv_params, "complex", ind = TRUE) # phase one data with biomarker indicator
      
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
      
      # influence function
      if_org = pkg_if(dat_org, t_org)
      if_ind = pkg_if(dat_ind, t_ind)
      
      # cox estimator with influence-function based inference
      val_n_tps_org = est_med(dat_org, t_org, edge = FALSE, boots = num_boot, if_result = if_org$tps)
      val_n_tps_ind = est_med(dat_ind, t_ind, edge = FALSE, boots = num_boot, if_result = if_ind$tps)
      val_n_flx_org = est_med(dat_org, t_org, edge = TRUE, boots = num_boot, if_result = if_org$flx)
      val_n_flx_ind = est_med(dat_ind, t_ind, edge = TRUE, boots = num_boot, if_result = if_ind$flx)
      
      # split results into X-only and X+S tables
      res_tps_org_raw = split_med_result(val_n_tps_org$result, "raw")
      res_tps_org_xs = split_med_result(val_n_tps_org$result, "xs")
      res_tps_ind_raw = split_med_result(val_n_tps_ind$result, "raw")
      res_tps_ind_xs = split_med_result(val_n_tps_ind$result, "xs")
      res_flx_org_raw = split_med_result(val_n_flx_org$result, "raw")
      res_flx_org_xs = split_med_result(val_n_flx_org$result, "xs")
      res_flx_ind_raw = split_med_result(val_n_flx_ind$result, "raw")
      res_flx_ind_xs = split_med_result(val_n_flx_ind$result, "xs")
      small_results = list(
        tps_org_raw = res_tps_org_raw,
        tps_org_xs = res_tps_org_xs,
        tps_ind_raw = res_tps_ind_raw,
        tps_ind_xs = res_tps_ind_xs,
        flx_org_raw = res_flx_org_raw,
        flx_org_xs = res_flx_org_xs,
        flx_ind_raw = res_flx_ind_raw,
        flx_ind_xs = res_flx_ind_xs
      )
      
      # results
      result = list(
        NIE_0_org = unname(val_0_org["NIE", "true"]),
        NIE_0_ind = unname(val_0_ind["NIE", "true"]),
        NDE_0_org = unname(val_0_org["NDE", "true"]),
        NDE_0_ind = unname(val_0_ind["NDE", "true"]),
        PM_0_org = unname(val_0_org["PM", "true"]),
        PM_0_ind = unname(val_0_ind["PM", "true"])
      )
      for (i in names(small_results)) {
        data_type = if (grepl("_org_", i)) {
          "org"
        } else {
          "ind"
        }
        truth = if (data_type == "org") {
          val_0_org
        } else {
          val_0_ind
        }
        result_table = small_results[[i]]
        
        for (j in c("NIE", "NDE", "PM")) {
          estimate = result_table["estimate", j]
          true_value = truth[j, "true"]
          # estimator
          result[[paste0(j, "_n_", i)]] = unname(estimate)
          # se from bootstrap or if
          result[[paste0(j, "_se_bs_", i)]] = unname(result_table["se_bs", j])
          result[[paste0(j, "_se_if_", i)]] = unname(result_table["se_if", j])
          # bias percentage
          result[[paste0(j, "_bias_pct_", i)]] = unname((estimate - true_value) / true_value * 100)
          # ci from bootstrap or if
          for (method in c("bs", "if")) {
            result[[paste0(j, "_low_", method, "_", i)]] = unname(result_table[paste0("low_", method), j])
            result[[paste0(j, "_up_", method, "_", i)]] = unname(result_table[paste0("up_", method), j])
          }
        }
      }
      result[[".complex"]] = list(
        dat_org = dat_org,
        dat_ind = dat_ind,
        val_0_org = val_0_org,
        val_0_ind = val_0_ind,
        val_n_tps_org = val_n_tps_org,
        val_n_tps_ind = val_n_tps_ind,
        val_n_flx_org = val_n_flx_org,
        val_n_flx_ind = val_n_flx_ind
      )
      
      return(result)
    })
  },
  
  main = {
    sim %<>% run()
    print(sim$errors)
  },
  
  last = {
    combos = c("tps_org_raw", "tps_org_xs", "tps_ind_raw", "tps_ind_xs", "flx_org_raw", "flx_org_xs", "flx_ind_raw", "flx_ind_xs")
    raw_combos = combos[grepl("_raw$", combos)]
    effects = c("NIE", "NDE", "PM")
    truth_name = function(i, j) {
      data_type = if (grepl("_org_", j)) {
        "org"
      } else {
        "ind"
      }
      paste0(i, "_0_", data_type)
    }
    summary_call = function(specs) {
      do.call(SimEngine::summarize, c(list(sim = sim), specs))
    }
    
    mean_specs = list()
    se_bs_specs = list()
    se_if_specs = list()
    bias_specs = list()
    bias_pct_specs = list()
    coverage_bs_specs = list()
    coverage_if_specs = list()
    
    for (j in combos) {
      for (i in effects) {
        estimate = paste0(i, "_n_", j)
        truth = truth_name(i, j)
        # mean of true values and estimators
        mean_specs[[length(mean_specs) + 1]] = list(stat = "mean", x = truth, name = paste0("mean_true_", i, "_", j))
        mean_specs[[length(mean_specs) + 1]] = list(stat = "mean", x = estimate, name = paste0("mean_est_", i, "_", j))
        # mean bootstrap se (all 8 combinations)
        se_bs_specs[[length(se_bs_specs) + 1]] = list(stat = "mean", x = paste0(i, "_se_bs_", j), name = paste0("se_bs_", i, "_", j))
        # bias
        bias_specs[[length(bias_specs) + 1]] = list(stat = "bias", estimate = estimate, truth = truth, name = paste0("bias_", i, "_", j))
        # bias percentage
        bias_pct_specs[[length(bias_pct_specs) + 1]] = list(stat = "mean", x = paste0(i, "_bias_pct_", j), name = paste0("bias_pct_", i, "_", j))
        # bootstrap coverage (all 8 combinations)
        coverage_bs_specs[[length(coverage_bs_specs) + 1]] = list(stat = "coverage", lower = paste0(i, "_low_bs_", j), upper = paste0(i, "_up_bs_", j), truth = truth, name = paste0("cov_bs_", i, "_", j))
        # influence function exists only for the 4 X-only combos
        if (j %in% raw_combos) {
          se_if_specs[[length(se_if_specs) + 1]] = list(stat = "mean", x = paste0(i, "_se_if_", j), name = paste0("se_if_", i, "_", j))
          coverage_if_specs[[length(coverage_if_specs) + 1]] = list(stat = "coverage", lower = paste0(i, "_low_if_", j), upper = paste0(i, "_up_if_", j), truth = truth, name = paste0("cov_if_", i, "_", j))
        }
      }
    }
    
    mean_results = summary_call(mean_specs)
    standard_error_bs = summary_call(se_bs_specs)
    standard_error_if = summary_call(se_if_specs)
    bias = summary_call(bias_specs)
    bias_percentage = summary_call(bias_pct_specs)
    coverage_bs = summary_call(coverage_bs_specs)
    coverage_if = summary_call(coverage_if_specs)
    
    summary_results = list(
      mean = mean_results,
      standard_error_bs = standard_error_bs,
      standard_error_if = standard_error_if,
      bias = bias,
      bias_percentage = bias_percentage,
      coverage_bs = coverage_bs,
      coverage_if = coverage_if
    )
  },
  
  cluster_config = list(js = "slurm")
)



# save results
# saveRDS(bias, file = "Evaluation/version_1/vaccine_bias.rds")
# saveRDS(bias_percentage, file = "Evaluation/version_1/vaccine_bias_percentage.rds")
# saveRDS(coverage, file = "Evaluation/version_1/vaccine_coverage.rds")
# saveRDS(estimators, file = "Evaluation/version_1/vaccine_estimaters.rds")
# saveRDS(standard_error, file = "Evaluation/version_1/vaccine_standard_error.rds")
# saveRDS(true_values, file = "Evaluation/version_1/vaccine_true_values.rds")



# end time
end_time = Sys.time()
execution_time = end_time - start_time
print(execution_time)





