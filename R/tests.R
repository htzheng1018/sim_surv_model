
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

#  calculate inverse-probability weights
dat_test$Ns = ave(rep(1, nrow(dat_test)), dat_test$trt, dat_test$cc_strata, FUN = length)
dat_test$ns = ave(dat_test$casecontrol, dat_test$trt, dat_test$cc_strata, FUN = sum)
dat_test$ipw = ifelse(dat_test$ns > 0, dat_test$Ns / dat_test$ns, 0)

# clean dataset
dat_test = dat_test %>%
  dplyr::transmute(Y = HIVwk28preunblfu,
                delta = HIVwk28preunbl,
                treat = trt,
                S1 = IgG_AEgp41ID_BP, S2 = IgG_BioRV144_C52_AE, S3 = IgG_BioRV144_C52_C, S4 = IgG_BioRV144_V2_AE, S5 = IgG_BioRV144_V2_B,
                X1 = age, X2 = BMI, X3 = bhvrisk,
                Z = casecontrol,
                ipw = ipw)
head(dat_test)



##############
#### test ####
##############

# for IgG_AEgp41ID_BP
dat_test1 = dat_test %>% dplyr::mutate(S = S1)
val_n_tps = est_med(dat_test1, 578)
val_n_flx = est_med(dat_test1, 578, edge = T)





