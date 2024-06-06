# we should use this: https://www.tidymodels.org/learn/statistics/survival-metrics/ 

library(survival)
library(tidyverse)
library(glmnet)

# Prep Data ---------------------------------------------------------------

load(here::here("inst", "cureit-simulation-results", "sim_data_list.RData"))
load(here::here("inst", "cureit-simulation-results", "sim_valid_data_list.RData"))
load(here::here("inst","cureit-simulation-results", "cv_fits_100.RData"))

dat <- sim_data_list[[1]]
fit <- cv_fits_list[[1]]
final_fit <- fit$fit[[fit$index$`min`[1]]][[fit$index$min[2]]]

get_brier_for_each_data <- function(dat, fit) {
    
  final_fit <- fit$fit[[fit$index$`min`[1]]][[fit$index$min[2]]]
  
  # Create Data frame to pass  -----
  
  data = data.frame(dat$x) %>%
    mutate(times = dat$t, 
           event = dat$d,
           "truth_survival" = Surv(dat$t, dat$d))
  
  # Predict Function (glmnet cure) ---------------------------------------------
  
  predict_df <- predict_cure(final_fit = final_fit, new_data = data)
  predict_df_long <- predict_df %>%
    unnest(everything())
  
  
  brier_scores_yardstick <-
    predict_df %>% 
    mutate(truth = Surv(times, event)) %>%
    brier_survival(truth = truth, .pred)
  
  return(brier_scores_yardstick)
  
}

# estimates greater than 1  
b <- get_brier_for_each_data(dat = sim_data_list[[3]], fit = cv_fits_list[[3]])


res <- list()
for (i in 1:10) {
  
  res[[i]] <- get_brier_for_each_data(sim_data_list[[i]], cv_fits_list[[i]])
}

ggplot(res[[3]], aes(x = .eval_time, .estimate)) + geom_point()




# quick brier ---------
dat <- sim_data_list[[3]]
fit <- cv_fits_list[[3]]
final_fit <- fit$fit[[fit$index$`min`[1]]][[fit$index$min[2]]]

# Create Data frame to pass  -----

data = data.frame(dat$x) %>%
  mutate(times = dat$t, 
         event = dat$d,
         "truth_survival" = Surv(dat$t, dat$d))

predict_df <- predict_cure(final_fit = final_fit, new_data = data)
predict_df_long <- predict_df %>%
  unnest(everything())


x <- c()
for (i in 1:100) {
  
  time_1 <- predict_df_long$.eval_time[[i]]
  
  predict_df_long_1 <- predict_df_long %>%
    filter(.eval_time == time_1)
  
  x[i] <-  mean(I(predict_df_long_1$times > time_1)*(1 - predict_df_long_1$.pred_survival)^2*(predict_df_long_1$.weight_censored) + 
                  I(predict_df_long_1$event == 1 & predict_df_long_1$times <= time_1)*
                  (0 - predict_df_long_1$.pred_survival)^2*(predict_df_long_1$.weight_censored))
  
  # OLD -----
  # x[i] <- mean(I(predict_df_long_1$event == 0)*(1 - predict_df_long_1$.pred_survival)^2/(predict_df_long_1$.weight_censored + 0.001) + 
  #                I(predict_df_long_1$event ==1)*(0 - predict_df_long_1$.pred_survival)^2/(predict_df_long_1$.weight_censored  + 0.001))
  # 
  # OLD ----
  # brier[l] <- mean(I(di==0)*(1 - preds)^2/(ipw+0.001) +
  #                    I(di==1)*(0 - preds)^2/(ipw+0.001))
  # 
  # NEW -----
  # brier[l] <- mean(I(ti > tbrier[l])*(1 - preds)^2/(ipw+0.001) + 
  #                    I(di==1 & ti <= tbrier[l])*(0 - preds)^2/(ipw+0.001))
  
}

times <- unique(predict_df_long$.eval_time)
brier <- x

brier_scores_manual <- bind_cols("times" = times, "brier" = brier)

brier_scores_yardstick <-
  predict_df %>% 
  mutate(truth = Surv(times, event)) %>%
  brier_survival(truth = truth, .pred)

# Brier With Yardstick -------------------------------------------------------
library(yardstick)


