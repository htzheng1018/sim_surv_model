
library(dplyr)
library(kableExtra)
summary_results = readRDS("Evaluation/summary_results.rds")

# save files
prepare_table = function(data, target) {
  data %>%
    select(n, surv_time, contains(target)) %>%
    select(-matches(paste0("NDE_", target, "_xs$"))) %>%
    mutate(surv_time = recode(surv_time, Exp = "Exponential", Gom = "Gompertz")) %>%
    rename(size = n, `survival type` = surv_time)
}
targets = c("tps_org", "tps_ind", "flx_org", "flx_ind")
for (target in targets) {
  prepare_table(summary_results$bias_percentage, target) %>%
    kable(format = "latex", booktabs = TRUE, longtable = FALSE, align = "r") %>%
    save_kable(paste0("/Users/gjp731/Desktop/research/Avi_Kenny/results/bias_percentage_", target, ".tex"))
  prepare_table(summary_results$coverage_bs, target) %>%
    mutate(across(starts_with("cov_"), ~ .x * 100)) %>%
    kable(format = "latex", booktabs = TRUE, longtable = FALSE, align = "r") %>%
    save_kable(paste0("/Users/gjp731/Desktop/research/Avi_Kenny/results/coverage_bs_", target, ".tex"))
  prepare_table(summary_results$coverage_if, target) %>%
    mutate(across(starts_with("cov_"), ~ .x * 100)) %>%
    kable(format = "latex", booktabs = TRUE, longtable = FALSE, align = "r") %>%
    save_kable(paste0("/Users/gjp731/Desktop/research/Avi_Kenny/results/coverage_if_", target, ".tex"))
  prepare_table(summary_results$standard_error_bs, target) %>%
    kable(format = "latex", booktabs = TRUE, longtable = FALSE, align = "r") %>%
    save_kable(paste0("/Users/gjp731/Desktop/research/Avi_Kenny/results/standard_error_bs_", target, ".tex"))
  prepare_table(summary_results$standard_error_if, target) %>%
    kable(format = "latex", booktabs = TRUE, longtable = FALSE, align = "r") %>%
    save_kable(paste0("/Users/gjp731/Desktop/research/Avi_Kenny/results/standard_error_if_", target, ".tex"))
}
