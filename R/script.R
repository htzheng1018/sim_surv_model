library(kableExtra)



# save files
kable(bias, format = "latex", booktabs = TRUE, longtable = FALSE) %>%
  save_kable("/Users/gjp731/Desktop/research/Avi_Kenny/results/bias.tex")
kable(bias_percentage, format = "latex", booktabs = TRUE, longtable = FALSE) %>%
  save_kable("/Users/gjp731/Desktop/research/Avi_Kenny/results/bias_percentage.tex")
kable(coverage, format = "latex", booktabs = TRUE, longtable = FALSE) %>%
  save_kable("/Users/gjp731/Desktop/research/Avi_Kenny/results/coverage.tex")
kable(estimators, format = "latex", booktabs = TRUE, longtable = FALSE) %>%
  save_kable("/Users/gjp731/Desktop/research/Avi_Kenny/results/estimators.tex")
kable(standard_error, format = "latex", booktabs = TRUE, longtable = FALSE) %>%
  save_kable("/Users/gjp731/Desktop/research/Avi_Kenny/results/standard_error.tex")
kable(true_values, format = "latex", booktabs = TRUE, longtable = FALSE) %>%
  save_kable("/Users/gjp731/Desktop/research/Avi_Kenny/results/true_values.tex")
