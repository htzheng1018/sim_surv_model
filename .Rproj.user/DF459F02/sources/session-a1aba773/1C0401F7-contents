
pkg_if = function(dat, t) {
  if (!requireNamespace("vaccine", quietly = TRUE)) {
    stop("Package 'vaccine' is required for influence-function inference.")
  }
  
  # Data for r_p and r_v_one: X is observed in the phase-one sample.
  dat_pkg = dat
  dat_pkg$S_pkg = ifelse(dat_pkg$treat == 1 & dat_pkg$Z == 0, NA_real_, dat_pkg$S)
  dat_pkg$S_pkg[dat_pkg$treat == 0] = 0
  dat_pkg$Z_pkg = ifelse(dat_pkg$treat == 0, 1, dat_pkg$Z)
  dat_pkg$ipw_pkg = ifelse(dat_pkg$treat == 0, 1, dat_pkg$ipw)
  
  dat_overall = vaccine::load_data(
    time = "Y",
    event = "delta",
    vacc = "treat",
    marker = "S_pkg",
    covariates = c("X1", "X2"),
    weights = "ipw_pkg",
    ph2 = "Z_pkg",
    data = dat_pkg,
    covariates_ph2 = FALSE
  )
  
  # Data for r_m: X and S enter through the phase-two/IPW estimator.
  dat_pkg$X1_pkg = ifelse(dat_pkg$treat == 1 & dat_pkg$Z == 0, NA_real_, dat_pkg$X1)
  dat_pkg$X2_pkg = ifelse(dat_pkg$treat == 1 & dat_pkg$Z == 0, NA_real_, dat_pkg$X2)
  
  dat_control = vaccine::load_data(
    time = "Y",
    event = "delta",
    vacc = "treat",
    marker = "S_pkg",
    covariates = c("X1_pkg", "X2_pkg"),
    weights = "ipw_pkg",
    ph2 = "Z_pkg",
    data = dat_pkg,
    covariates_ph2 = TRUE
  )
  
  # vaccine only returns these vectors when this internal option is TRUE.
  vaccine_env = getFromNamespace(".vaccine_env", "vaccine")
  old_return_IF_vec = vaccine_env$return_IF_vec
  on.exit(
    vaccine_env$return_IF_vec <- old_return_IF_vec,
    add = TRUE
  )
  vaccine_env$return_IF_vec = TRUE
  
  # This part is shared by the TPS and FLX models.
  overall = vaccine::est_overall(
    dat = dat_overall,
    t_0 = t,
    method = "Cox",
    ve = FALSE
  )
  
  r_p = overall$est[overall$group == "placebo"]
  r_v = overall$est[overall$group == "vaccine"]
  
  # est_overall stores the IF of survival. Risk equals one minus survival.
  IF_p = -attr(overall, "IF_vec_placebo")
  IF_v = -attr(overall, "IF_vec_vaccine")
  
  n_p = sum(dat$treat == 0)
  n_v = sum(dat$treat == 1)
  
  get_if = function(edge) {
    controlled = vaccine::est_ce(
      dat = dat_control,
      type = "Cox",
      t_0 = t,
      cr = TRUE,
      cve = FALSE,
      s_out = 0,
      ci_type = "regular",
      params_cox = vaccine::params_ce_cox(edge_ind = edge)
    )
    
    r_m = controlled$cr$est[1]
    IF_m = attr(controlled, "IF_vec_rM")
    
    if (length(IF_v) != n_v ||
        length(IF_m) != n_v ||
        length(IF_p) != n_p) {
      stop("Unexpected influence-vector lengths returned by vaccine.")
    }
    
    NIE = r_v / r_m
    NDE = r_m / r_p
    PM = 1 - log(NDE) / (log(NIE) + log(NDE))
    
    # Delta method for NIE = r_v / r_m.
    IF_NIE = IF_v / r_m - r_v * IF_m / r_m ^ 2
    var_NIE = sum(IF_NIE ^ 2) / n_v ^ 2
    
    # Delta method for NDE = r_m / r_p.
    IF_NDE_v = IF_m / r_p
    IF_NDE_p = -r_m * IF_p / r_p ^ 2
    
    var_NDE = sum(IF_NDE_v ^ 2) / n_v ^ 2 +
      sum(IF_NDE_p ^ 2) / n_p ^ 2
    
    # Delta method for
    # PM = 1 - log(r_m / r_p) / log(r_v / r_p).
    A = log(r_m / r_p)
    B = log(r_v / r_p)
    
    IF_PM_v = A * IF_v / (B ^ 2 * r_v) -
      IF_m / (B * r_m)
    IF_PM_p = log(r_v / r_m) * IF_p / (B ^ 2 * r_p)
    
    var_PM = sum(IF_PM_v ^ 2) / n_v ^ 2 +
      sum(IF_PM_p ^ 2) / n_p ^ 2
    
    estimate = c(
      NIE_one = NIE,
      NDE = NDE,
      PM_one = PM
    )
    
    se = sqrt(c(
      NIE_one = var_NIE,
      NDE = var_NDE,
      PM_one = var_PM
    ))
    
    ratio_ci = function(est, se) {
      if (!is.finite(est) || !is.finite(se) || est <= 0) {
        return(c(NA_real_, NA_real_))
      }
      
      exp(log(est) + c(-1, 1) * 1.96 * se / est)
    }
    
    ci_NIE = ratio_ci(NIE, se["NIE_one"])
    ci_NDE = ratio_ci(NDE, se["NDE"])
    ci_PM = PM + c(-1, 1) * 1.96 * se["PM_one"]
    
    result = rbind(
      estimate = estimate,
      se_if = se,
      low_if = c(
        NIE_one = unname(ci_NIE[1]),
        NDE = unname(ci_NDE[1]),
        PM_one = unname(ci_PM[1])
      ),
      up_if = c(
        NIE_one = unname(ci_NIE[2]),
        NDE = unname(ci_NDE[2]),
        PM_one = unname(ci_PM[2])
      )
    )
    
    return(result)
  }
  
  # r_p and r_v are computed once and reused here.
  return(list(
    tps = get_if(FALSE),
    flx = get_if(TRUE)
  ))
}

est_med = function(dat, t, edge = FALSE, boots = 1000, if_result = NULL) {
  dat$S_zero = as.numeric(dat$S == 0) # indicator of S = 0
  X_all = grep("^X", names(dat), value = TRUE)
  
  # survival function formular
  form_plc = as.formula(paste("Surv(Y, delta) ~", paste(X_all, collapse = " + ")))
  if (edge == FALSE) {
    form_vac_two = as.formula(paste("Surv(Y, delta) ~", paste(X_all, collapse = " + "), "+ S"))
  } else {
    form_vac_two = as.formula(paste("Surv(Y, delta) ~", paste(X_all, collapse = " + "), "+ S + S_zero"))
  }
  
  # estimated risk ratio
  risk_n = function(model, data, treatment) {
    bh = basehaz(model, centered = FALSE)
    # cumulative baseline hazard
    before_t = which(bh$time <= t)
    H_0 = if (length(before_t) == 0) {
      0
    } else {
      bh$hazard[max(before_t)]
    }
    # estimated survival functions
    beta = model$coefficients
    X_S = data[, names(beta), drop = FALSE]
    unprop = drop(exp(as.matrix(X_S) %*% beta))
    Q = exp(-H_0 * unprop)
    if (treatment %in% c("plc", "vac_one")) {
      Q_n = mean(Q)
    } else if (treatment == "vac_two") {
      Q_n = weighted.mean(Q, w = data$ipw)
    } else if (treatment == "med") {
      X_S_med = X_S
      X_S_med$S = 0
      if (edge == TRUE) {X_S_med$S_zero = 1} # indicator of S= 0
      unprop_med = drop(exp(as.matrix(X_S_med) %*% beta))
      Q_med = exp(-H_0 * unprop_med)
      Q_n = weighted.mean(Q_med, w = data$ipw)
    }
    return(1 - Q_n)
  }
  
  get_estimates = function(dat) {
    # get data in different treatment groups
    dat_plc = dat[dat$treat == 0, , drop = FALSE]
    dat_vac_one = dat[dat$treat == 1, , drop = FALSE]
    dat_vac_two = dat[dat$treat == 1 & dat$Z == 1, , drop = FALSE]
    
    # coxph models
    model_plc = coxph(form_plc, data = dat_plc, model = TRUE)
    model_vac_one = coxph(form_plc, data = dat_vac_one, model = TRUE)
    model_vac_two = coxph(form_vac_two, data = dat_vac_two, weights = ipw, model = TRUE)
    
    # estimaton of NIE, NDE and proportion mediated
    r_p = risk_n(model_plc, dat_plc, "plc")
    r_v_one = risk_n(model_vac_one, dat_vac_one, "vac_one")
    r_v_two = risk_n(model_vac_two, dat_vac_two, "vac_two")
    r_m = risk_n(model_vac_two, dat_vac_two, "med")
    NIE_one = r_v_one / r_m
    NIE_two = r_v_two / r_m
    NDE = r_m / r_p
    PM_one = 1 - log(NDE) / (log(NIE_one) + log(NDE))
    PM_two = 1 - log(NDE) / (log(NIE_two) + log(NDE))
    
    return(list(
      estimate = c(NIE_one = NIE_one, NIE_two = NIE_two, NDE = NDE, PM_one = PM_one, PM_two = PM_two),
      risks = c(r_p = r_p, r_v_one = r_v_one, r_v_two = r_v_two, r_m = r_m)
    ))
  }
  
  # results of estimates
  fit = get_estimates(dat)
  result = matrix(NA_real_, nrow = 7, ncol = 5, dimnames = list(
    c("estimate", "se_bs", "low_bs", "up_bs", "se_if", "low_if", "up_if"),
    c("NIE_one", "NIE_two", "NDE", "PM_one", "PM_two")
  ))
  result["estimate", ] = fit$estimate
  
  # bootstrap for ci
  if (boots > 0) {
    boot_est = matrix(NA_real_, nrow = boots, ncol = length(fit$estimate),
                      dimnames = list(NULL, names(fit$estimate))
    )
    for (i in seq_len(boots)) {
      samps = sample.int(nrow(dat), size = nrow(dat), replace = TRUE)
      boot_est[i, ] = tryCatch(get_estimates(dat[samps, , drop = FALSE])$estimate,
                               error = function(e) {
                                 rep(NA_real_, length(fit$estimate))
                               }
      )
    }
    n_boot_valid = colSums(is.finite(boot_est))
    boot_se = vapply(seq_len(ncol(boot_est)),
                     function(j) {
                       x = boot_est[, j]
                       x = x[is.finite(x)]
                       if (length(x) < 2) {
                         NA_real_
                       } else {
                         sd(x)
                       }
                     },
                     numeric(1)
    )
    boot_ci = vapply(seq_len(ncol(boot_est)),
                     function(j) {
                       x = boot_est[, j]
                       x = x[is.finite(x)]
                       if (length(x) == 0) {
                         c(NA_real_, NA_real_)
                       } else {
                         quantile(x, probs = c(0.025, 0.975), names = FALSE)
                       }
                     },
                     numeric(2)
    )
    result["se_bs", ] = boot_se
    result["low_bs", ] = boot_ci[1, ]
    result["up_bs", ] = boot_ci[2, ]
  }
  
  # influence function for ci
  if (is.null(if_result)) {
    all_if = pkg_if(dat, t)
    if_result = if (edge == FALSE) {
      all_if$tps
    } else {
      all_if$flx
    }
  }
  if_names = c("NIE_one", "NDE", "PM_one")
  result[c("se_if", "low_if", "up_if"), if_names] = if_result[c("se_if", "low_if", "up_if"), if_names, drop = FALSE]
  
  # results of risks
  risks = data.frame(t(fit$risks), check.names = FALSE)
  
  return(list(result = result, risks = risks))
}
