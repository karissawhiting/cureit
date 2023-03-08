library(tidycmprsk)

x <- tidycmprsk::trial %>% select(trt,
                      age,
                      ttdeath,
                      death_cr
) %>%
  tbl_uvregression(method = tidycmprsk::crr,
                   y = Surv(ttdeath, death_cr),
                   exponentiate = TRUE) %>%
                   method.args = list(failcode = c("death from cancer",
                                                   "death other causes")) %>%
  bold_labels() %>%
  bold_p()
x




m1 <- glmmTMB::glmmTMB(count ~ mined,
                       ziformula = ~ mined,
                       family = poisson,
                       data = glmmTMB::Salamanders)
broom.mixed::tidy_plus_plus(m1)
tbl_regression(m1)

tbl_uvregression()

multi2  <- cureit(formula = Surv(os_months, dod) ~ sex_dscrp,
                  cure_formula = ~ sex_dscrp,
                  data = df_cure,
                  nboot = 0)

multi2 %>%
  broom.helpers::tidy_plus_plus()
broom.mixed::tidy(multi2)




x <- cureit(formula = Surv(ttdeath, death) ~ age + grade,
            cure_formula = ~ age + grade,  data = trial)



y <- cureit(formula = Surv(ttdeath, death) ~ age + grade,
  cure_formula = NULL,  data = trial)

identical(x, y)
