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

# Create Data frame to pass 

data = data.frame(dat$x) %>%
  mutate(times = dat$t, 
         event = dat$d,
         "truth_survival" = Surv(dat$t, dat$d))

# * Calculate Sensitivity of Variable Selection in Training Fits -------
spe_cox_min <- c()
sen_cox_min <- c()

x <- for (i in 1:100) {
  fit <- cv_fits_list[[i]]
  #fit <- x$fit
  dat <- sim_data_list[[i]]
  
  spe_cox_min[[i]] <- 1-mean(setdiff(1:200,dat$id_cure) %in% which(fit$fit[[fit$index$min[1]]][[fit$index$min[2]]]$beta!=0))
  sen_cox_min[[i]] <- mean(dat$id_cox %in% which(fit$fit[[fit$index$min[1]]][[fit$index$min[2]]]$beta!=0))
  # spe_cox_min <- 1-mean(setdiff(1:200,dat$id_cure) %in% which(fit$fit[[fit$index$min[1]]][[fit$index$min[2]]]$beta!=0))
  
  # print(sen_cox_min)
  
}


hist(unlist(sen_cox_min))
hist(unlist(spe_cox_min))

# Check Variable Selection on Standard Lasso
# will add foldid (see new version of cv.lassocure())
cv.glmnet(x = dat$x, y = Surv(dat$t, dat$d), family = "cox")

# Predict Function (glmnet cure) ---------------------------------------------

# Trying to get predictions in similar format as tidymodels yardstick
# Eventually will integrate this with existing predict method function



predict_df <- predict_cure(final_fit = final_fit, new_data = data)
predict_df_long <- predict_df %>%
  unnest(everything())



# quick brier ---------
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

# Brier With Yardstick -------------------------------------------------------
library(yardstick)

# Doesn't match- weights are wrong
brier_scores_yardstick <-
  predict_df %>% 
  mutate(truth = Surv(times, event)) %>%
  brier_survival(truth = truth, .pred)


all_brier <- brier_scores_yardstick %>%
  select(times = .eval_time, brier_yardstick = .estimate) %>%
  left_join(brier_scores_manual)

# Brier Function -------------------------------------------------------

predict_df_long <- predict_df %>%
  unnest(everything())

predict_df_nest_time <- predict_df_long %>%
  nest(data = -.eval_time)

predict_df_timepoint = predict_df_nest_time[10,] %>%
  unnest(everything())

brier_for_each_eval <- function(predict_df_timepoint) {
  
  
  x = mean(I(predict_df_timepoint$times > predict_df_timepoint$.eval_time )*(1 - predict_df_timepoint$.pred_survival)^2/(predict_df_timepoint$.weight_censored+0.001) + 
                  I(predict_df_timepoint$event == 1 & predict_df_timepoint$times <= predict_df_timepoint$.eval_time  )*
                  (0 - predict_df_timepoint$.pred_survival)^2/(predict_df_timepoint$.weight_censored+0.001))
  x
  }
  

brier_for_each_eval <- function(predict_df_timepoint) {
  
  
  x = mean(
    I(predict_df_timepoint$times > predict_df_timepoint$.eval_time)*(1 - predict_df_timepoint$.pred_survival)^2/(predict_df_timepoint$.weight_censored+0.001) + 
      
             I(predict_df_timepoint$event == 1 & predict_df_timepoint$times <= predict_df_timepoint$.eval_time  )*
             (0 - predict_df_timepoint$.pred_survival)^2/(predict_df_timepoint$.weight_censored+0.001)
    )
  x/100
}

brier_for_each_eval2 <- function(predict_df_timepoint) {
  
  category1 = I(predict_df_timepoint$event == 1 & predict_df_timepoint$times <= predict_df_timepoint$.eval_time  )
  category2 = I(predict_df_timepoint$times > predict_df_timepoint$.eval_time)
  
  x = mean(
    category1*(0 - predict_df_timepoint$.pred_survival)^2/(predict_df_timepoint$.weight_censored) +
    category2 *(1 - predict_df_timepoint$.pred_survival)^2/(predict_df_timepoint$.weight_censored) 
     
  )
  x

}

brier_for_each_eval(predict_df_timepoint)

# Compare to yardstick ----

brier_survival_impl <- function(truth,
                                estimate,
                                censoring_weights,
                                case_weights,
                                eval_time) {
  surv_time <- yardstick:::.extract_surv_time(truth)
  surv_status <- yardstick:::.extract_surv_status(truth)
  
  if (!is.null(case_weights)) {
    case_weights <- vec_cast(case_weights, to = double())
    norm_const <- sum(case_weights)
  } else {
    case_weights <- rep(1, length(estimate))
    norm_const <- sum(!survival::is.na.Surv(truth))
  }
  
  category_1 <- surv_time <= eval_time & surv_status == 1
  category_2 <- surv_time > eval_time
  
  # (0 - estimate) ^ 2 == estimate ^ 2
  res <- (category_1 * estimate^2 * censoring_weights) +
    (category_2 * (1 - estimate)^2 * censoring_weights)
  
  res <- res * case_weights
  res <- sum(res, na.rm = TRUE)
  res / norm_const
}
predict_df_timepoint2 <- predict_df_timepoint %>%
  mutate(surv_obj = Surv(predict_df_timepoint$times, predict_df_timepoint$event))

library(yardstick)
brier_survival_impl(
                    truth = predict_df_timepoint2$surv_obj,
                    
                    estimate = predict_df_timepoint2$.pred_survival,
                                censoring_weights = predict_df_timepoint2$.weight_censored,
                                case_weights = NULL,
                                eval_time = predict_df_timepoint2$.eval_time)

brier_for_each_eval(predict_df_timepoint)

# Simple Former Brier Calculation ---------------------------------------------------
library(glmnet)
library(tidyverse)
library(survival)

calc_brier_old <- function(new_data, fit) {
  # ti <- dat_valid$t
  # di <- dat_valid$d
  # xi <- dat_valid$x
  ti <- new_data$t
  di <- new_data$d
  xi <- new_data$x
  
  # tj from train data
  tj <- fit$fit[[fit$index$`min`[1]]][[fit$index$`min`[2]]]$tj
  final_fit <- fit$fit[[fit$index$`min`[1]]][[fit$index$`min`[2]]]
  
  fitcure <- final_fit$fitcure
  predcure <- predict(fitcure,newx=xi,s=min(fitcure$lambda),type="response")
  
  fitcox <- final_fit$fitcox
  predsurvexp <- predict(fitcox,newx=xi,s=min(fitcox$lambda),type="response")
  
  haz <- final_fit$haz
  cumhaz <- final_fit$cumhaz
  
  # censoring weights from test data
  sf <- survfit(Surv(ti,1-di)~1)
  sf_pred <- summary(sf, time = ti)
  
  predcens <- sf_pred$surv
  tcens <- sf_pred$time
  
  # tbrier is evaluation timepoints (default to all timepoints in test data)
  tbrier <- sort(ti)
  brier <- rep(NA,length(tbrier))
  ipw <- rep(NA,length(tbrier))
  
  for (l in 1:length(tbrier)){
    
    predsurv <- rep(NA,length(ti)) # Conditional survival prob for all patients at tbrier[l]
    
    if (tbrier[l] >= min(tj)){
      
      ids <- which(tj == max(tj[tj <= tbrier[l]]))
      predsurv <- exp(-cumhaz[ids]*predsurvexp)
      ipw <- as.numeric(ti > tbrier[l])*predcens[l] + as.numeric(ti <= tbrier[l])*predcens
      
    }else if (tbrier[l] < min(tj)){
      
      predsurv <- rep(1,length(ti))
      ipw <- 1
      
    }
    
    preds <- 1 - predcure + predcure*predsurv
    brier[l] <- mean(I(di==0)*(1 - preds)^2/(ipw+0.001) + I(di==1)*(0 - preds)^2/(ipw+0.001))
    #ipw[l] <- ipw
   # print(ipw)
  }
  
  preds <- cbind.data.frame(tbrier, brier)
  return(preds)
 # plot(tbrier,brier,type="S")
  
}

x <- calc_brier_old(sim_data_list[[1]], cv_fits_list[[1]])
x_val <- calc_brier_old(sim_valid_data_list[[1]], cv_fits_list[[1]])

x
plot(x$tbrier,x$brier,type="S")
plot(x_val$tbrier,x_val$brier,type="S")


all_brier <- list()

for (i in 1:length(sim_data_list)) {
  all_brier[[i]] <- calc_brier_old( new_data = sim_data_list[[i]], fit = cv_fits_list[[i]])
}

all_brier_valid <- list()

for (i in 1:length(sim_valid_data_list)) {
  all_brier_valid[[i]] <- calc_brier_old( new_data = sim_valid_data_list[[i]],
                                          fit = cv_fits_list[[i]])
}

x <- all_brier_valid[[1]]
x_valid <- all_brier[[1]]

ggplot(x, aes(x =  tbrier,y = brier)) + geom_line( aes(x =  tbrier,y = brier)) +
  geom_line(data = x_valid, aes(x =tbrier, y = brier), color = "blue")

plot(x$tbrier,x$brier,type="S")
line(x_valid$tbrier,x_valid$brier,type="S")


traps <- c()

for (i in 1:length(all_brier)) {
  traps[i] <- trapz(all_brier[[i]]$tbrier,all_brier[[i]]$brier)

}


all_brier2 <- all_brier %>%
  do.call("rbind", .)

x <- all_brier2 %>% group_by(round(tbrier, 1)) %>%
  summarize(mean = mean(brier))

plot(x$`round(tbrier, 1)`, x$mean,type="S")



