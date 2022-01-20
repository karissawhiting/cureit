# library(smcure)
# library(ISwR)
# library(gtsummary)
# library(tidyverse)
# 
# # Prepare Data --------------------------
# mel <- ISwR::melanom %>%
#   mutate(status = case_when(status %in% c(1, 3) ~ 1,
#                             TRUE ~ 0)) %>%
#   mutate(thick_cat =
#               case_when(thick <= 129 ~ "≤ 129",
#                         thick <= 322 ~ "less equal 322",
#                         thick > 322 ~ "greater 322")) %>%
#   mutate(thick_cat = fct_relevel(thick_cat, "≤ 129", "≤ 322")) %>%
#   mutate(sex = case_when(sex == 2 ~ "male",
#                          sex == 1 ~ "female"),
#          ulc = case_when(ulc == 1 ~ "present",
#                          ulc == 2 ~ "absent"))
# 
# 
# x <- prepare_variables(formula = Surv(days, status) ~ ulc + sex + thick_cat + thick,
#                   data = mel)
# 
# x <- prepare_variables(formula = Surv(days, status) ~ thick,
#                   data = mel)
# 
# data(e1684)
# # fit PH mixture cure model
# pd <- smcure(Surv(FAILTIME,FAILCENS)~TRT+SEX+AGE,cureform=~TRT+SEX+AGE,
#      data=e1684,
#      model="ph",Var = FALSE)
# 
# 
# x <- prepare_variables(formula = Surv(FAILTIME,FAILCENS)~TRT+SEX+AGE,
#                   data = e1684)
# 
# # Try Parallel --------------------------
# 
# d <- e1684[1:150,]
# 
# library(tictoc)
# 
# #1170.896/60 =19.51493
# tic()
# test2 <- fit_cure(formula = Surv(FAILTIME,FAILCENS)~TRT + AGE,
#                   data = d,
#                   eps = .01, nboot = 20, n_runs = 3)
# toc()
# 
# # 520.561 /60
# tic()
# test3 <- fit_cure(formula = Surv(FAILTIME,FAILCENS)~TRT + AGE,
#                   data = d,
#                   eps = .01, nboot = 20, n_runs = 3,
#                   parallel = TRUE, workers = 3)
# toc()
# 
# 
# test4 <- fit_cure(formula = Surv(FAILTIME,FAILCENS)~TRT + AGE,
#                   data = d,
#                   eps = .01, nboot = 20, n_runs = 1,
#                   parallel = TRUE, workers = 3)
# toc()
# 
# cure_nomogram(fit, prediction_time = 1000)
# 
# cure_calibration(fit, prediction_time = 1000)
# 
# # Clean Multivariate
# 
# load("/Users/Whiting/Repositories/msk-project-repos/canter-2008-cure-models/outputs/fitted-models/fit_multi.RData")
# 
# fit_multi %>%
#   tidy_cure()
# 
# 
# 
# fit_multi %>%
#   tidy_cure_frac() %>% View()
# 
# 
# fit_multi %>%
#   tidy_surv() %>% View()
# 



# 
# tic()
# fit_multi <- fit_cure(formula = Surv(tt.DOD, DOD) ~ primary_site + primary_size_category + r_status, 
#          data = data, 
#          nboot = 50,
#          eps = .01,
#          n_runs = 3,
#          parallel = TRUE,
#          workers = 3)
# toc()
# beep(5)
