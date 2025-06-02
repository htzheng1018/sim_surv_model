library(kableExtra)



# save files
kable(bias_Q_pct_two, format = "latex", booktabs = TRUE, longtable = FALSE) %>%
  save_kable("/Users/haotianzheng/Desktop/research/Avi_Kenny/results/bias_Q_pct_two.tex")
kable(bias_Q_pct_km, format = "latex", booktabs = TRUE, longtable = FALSE) %>%
  save_kable("/Users/haotianzheng/Desktop/research/Avi_Kenny/results/bias_Q_pct_km.tex")

kable(coverage_two, format = "latex", booktabs = TRUE, longtable = FALSE) %>%
  save_kable("/Users/haotianzheng/Desktop/research/Avi_Kenny/results/coverage_two.tex")
kable(coverage_km, format = "latex", booktabs = TRUE, longtable = FALSE) %>%
  save_kable("/Users/haotianzheng/Desktop/research/Avi_Kenny/results/coverage_km.tex")