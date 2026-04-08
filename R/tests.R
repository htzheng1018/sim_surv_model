
################
#### set up ####
################

library(SimEngine)
library(kableExtra)
library(readr)
source("R/create_data.R", local = T)
source("R/true_func.r", local = T)
source("R/est_med.r", local = T)



###################
#### load data ####
###################

dat_primary = read.csv("data/primary505_for_sharing_upd.csv", stringsAsFactors = FALSE)
dat_bama = read.csv("data/bama.m.for_sharing.csv", stringsAsFactors = FALSE)

vrb_interest = c("pub_id", "IgG_AEgp41ID_BP", "IgG_BioRV144_C52_AE", "IgG_BioRV144_C52_C", "IgG_BioRV144_V2_AE", "IgG_BioRV144_V2_B")
dat_test = dat_bama[, vrb_interest, drop = FALSE]
dat_test = merge(dat_primary, dat_test, by = "pub_id", all.x = TRUE)
dat_test = dat_test[which(dat_test$week28 == 1 & !is.na(dat_test$BMI)), ]

print(paste("dimension of the dataset:", nrow(dat_test))) # 2303



##############
#### test ####
##############







