library(smcure)
library(ISwR)
library(gtsummary)
library(tidyverse)

# Prepare Data --------------------------
mel <- ISwR::melanom %>%
  mutate(status = case_when(status %in% c(1, 3) ~ 1, 
                            TRUE ~ 0)) %>%
  mutate(thick_cat = 
              case_when(thick <= 129 ~ "≤ 129", 
                        thick <= 322 ~ "less equal 322", 
                        thick > 322 ~ "greater 322")) %>%
  mutate(thick_cat = fct_relevel(thick_cat, "≤ 129", "≤ 322")) %>%
  mutate(sex = case_when(sex == 2 ~ "male", 
                         sex == 1 ~ "female"),
         ulc = case_when(ulc == 1 ~ "present", 
                         ulc == 2 ~ "absent"))

# Run Pipeline --------------------------
fit <- fit_cure(formula = Surv(days, status) ~ ulc + sex + thick_cat + thick + no,
                data = mel)

cure_nomogram(fit, prediction_time = 300)

cure_calibration(fit, prediction_time = 300)

stability <- multiple_mod_runs(formula = Surv(days, status) ~ ulc + sex + thick_cat + thick + no,
                data = mel)

stability
