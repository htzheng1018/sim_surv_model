
surv_type = "Gompertz"
surv_params = c(0.1, 1e-3)

dat_org = create_data(1000, surv_type, surv_params, "complex")
t_org = 40
val_0_org = true_func(surv_type, surv_params, t_org, dat_org, "math")
val_0_org
val_n_tps_org = est_med_if(dat_org, t_org, edge = FALSE)
val_n_tps_org
val_local_tps_org = est_med(dat_org, t_org, boots = 0)$result
val_local_tps_org
