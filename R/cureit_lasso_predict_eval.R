# we should use this: https://www.tidymodels.org/learn/statistics/survival-metrics/ 

library(survival)
library(tidyverse)
library(glmnet)

# Prep Data ---------------------------------------------------------------

load(here::here("inst", "cureit-simulation-results", "sim_data_list.RData"))
load(here::here("inst", "cureit-simulation-results", "sim_valid_data_list.RData"))
dat <- sim_data_list[[1]]

load(here::here("inst","cureit-simulation-results", "cv_fits_100.RData"))
fit <- cv_fits_list[[1]]
final_fit <- fit$fit[[fit$index$`min`[1]]][[fit$index$min[2]]]

# Create Data frame to pass 

data = data.frame(dat$x) %>%
  mutate(times = dat$t, 
         event = dat$d,
         "truth_survival" = Surv(dat$t, dat$d))

# Predict Function (glmnet cure) ---------------------------------------------

# Trying to get predictions in similar format as tidymodels yardstick
# Eventually will integrate this with existing predict method function

predict_cure <- function(final_fit,
                         data = NULL,
                         eval_timepoints = NULL) {

  # Get values from observed data
  survival_times = data$times
  sorted_survival_times <- sort(survival_times)
  sorted_event_times = sort(data$times[data$event==1])
  censoring_weights = survfit(Surv(data$times, 1 - data$event) ~ 1)$surv
  
  data_x <- data %>%
    select(-c(times, event, truth_survival)) %>%
    as.matrix()
  
  
  # Extract cumulative hazard from final fit (for event times)
  haz <- final_fit$haz
  cumhaz <- final_fit$cumhaz
  
  # default eval timepoints are observed event times if not supplied by user
  if(is.null(eval_timepoints)) {
    eval_timepoints = sorted_survival_times
  }


  # Get Relative Cure and Cox Risk (independent of timepoint) -------------------

  # Cure probability risk prediction
  fitcure <- final_fit$fitcure
  predcure <- predict(fitcure, newx = data_x,
                      s = min(fitcure$lambda), type = "response")
  
  # Cox survival probability risk prediction
  fitcox <- final_fit$fitcox
  predsurvexp <- predict(fitcox, 
                         newx = data_x,
                         s = min(fitcox$lambda),
                         type = "response")


  
  # Get Predicted probabilities at Selected Timepoints -------------------
  
  all_preds_all_tps <- list()
  
  # for each eval timepoint, get predicted survival values
  for (l in 1:length(eval_timepoints)) {
    all_preds <- data.frame("eval_timepoint" = NA, 
                            "preds" = NA, 
                            "pred_cure" = NA, 
                            "preds_surv" = NA, 
                            "weights" = NA)
    
    # vector for conditional survival prob for all patients at sorted_survival_times[l]
    #predsurv <- rep(NA,length(survival_times)) 
    

    if (eval_timepoints[l] >= min(sorted_event_times)){
      
      # maximum event time that is less than timepoint we are evaluating 
      ids <- which(
        sorted_event_times == 
          max(sorted_event_times[sorted_event_times <= eval_timepoints[l]])
        )
      
      # Conditional survival prob for all patients at given time
      predsurv <- exp(-cumhaz[ids]*predsurvexp)
        
      # inverse probability of the censoring weights
      ipw <- as.numeric(survival_times > eval_timepoints[l])*censoring_weights[l] + 
        as.numeric(survival_times <= eval_timepoints[l])*censoring_weights
      
      # if no events occurred before eval time 
    } else if (eval_timepoints[l] < min(sorted_event_times)) {
      
      predsurv <- rep(1,length(survival_times))
      ipw <- 1
      
    }
    
    preds <- 1 - predcure + predcure*predsurv
    
    all_preds <- tibble(
      
      ".eval_time" = eval_timepoints[l],
      ".pred_survival" = preds, 
      ".predcure" = predcure, 
      ".predsurv" = predsurv, 
      ".weight_censored_raw" = censoring_weights[l],
      ".weight_censored" = ipw) %>%
      mutate(id = 1:nrow(.))
    # 
    # names(all_preds) <- c(
    #                       "preds",
    #                       "predcure", "predsurv", "weights",
    #                       "id_num")
    
    all_preds_nest <-  all_preds %>% nest(.preds = -c(id))
    
    all_preds_all_tps[[l]] <- all_preds
                            
  #  names(all_preds) = c("pred", "pred_cure","pred_surv")
  }
  
  
  all_pred_df <- all_preds_all_tps %>% 
    do.call("rbind", .)
  
  all_pred_df <- all_pred_df %>%
    nest(.pred = -c(id)) %>% 
    bind_cols(., "times" = data$times, 
              "event" = data$event)

  
  return(all_pred_df)
}

predict_df <- predict_cure(final_fit = final_fit, data = data)

# quick brier
x <- c()
for (i in 1:10) {
  
  time_1 <- predict_df_long$.eval_time[[i]]
  
  predict_df_long_1 <- predict_df_long %>%
    filter(.eval_time == time_1)
  
  x[i] <- mean(I(predict_df_long_1$event == 0)*(1 - predict_df_long_1$.pred_survival)^2/(predict_df_long_1$.weight_censored + 0.001) + 
         I(predict_df_long_1$event ==1)*(0 - predict_df_long_1$.pred_survival)^2/(predict_df_long_1$.weight_censored  + 0.001))
}

# Brier With Yardstick -------------------------------------------------------
library(yardstick)

brier_scores <-
  predict_df %>% 
  mutate(truth = Surv(times, event)) %>%
  brier_survival(truth = truth, .pred)


# Brier Function -------------------------------------------------------

brier_cure <- function(predict_df, truth = times, eval_timepoints = NULL) {

  
  # Get values from observed data
  survival_times = predict_df$times
  sorted_survival_times <- sort(survival_times)
  sorted_event_times = sort(data$times[data$event==1])
  
  
  predict_df_long <- predict_df %>%
    unnest(everything())
  
  predict_df_nest_time <- predict_df_long %>%
    nest(data = -.eval_time)
  

  map(predict_df_nest_time$data, function(x) {
    mean(I(x$event == 0)*(1 - x$.pred_survival)^2/(x$.weight_censored + 0.001) + I(x$event == 1)*(0 - x$.pred_survival)^2/(x$.weight_censored + 0.001))
    
  })
  for (l in 1:length(eval_timepoints)) {
    brier[l] <- mean(I(di==0)*(1 - preds)^2/(ipw+0.001) + I(di==1)*(0 - preds)^2/(ipw+0.001))
  
  }
}
  


# Former Brier Calculation ---------------------------------------------------

calc_brier_old <- function(dat, fit) {
  # ti <- dat_valid$t
  # di <- dat_valid$d
  # xi <- dat_valid$x
  ti <- dat$t
  di <- dat$d
  xi <- dat$x
  tj <- fit$fit[[fit$index$`min`[1]]][[fit$index$`min`[2]]]$tj
  final_fit <- fit$fit[[fit$index$`min`[1]]][[fit$index$`min`[2]]]
  
  fitcure <- final_fit$fitcure
  predcure <- predict(fitcure,newx=xi,s=min(fitcure$lambda),type="response")
  fitcox <- final_fit$fitcox
  predsurvexp <- predict(fitcox,newx=xi,s=min(fitcox$lambda),type="response")
  haz <- final_fit$haz
  cumhaz <- final_fit$cumhaz
  
  predcens <- survfit(Surv(ti,1-di)~1)$surv
  tcens <- survfit(Surv(ti,1-di)~1)$time
  
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

for (i in 1:length(cv_fits_list)) {
  
}
