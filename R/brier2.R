# we should use this: https://www.tidymodels.org/learn/statistics/survival-metrics/ 

library(survival)
library(tidyverse)

# Prep Data ---------------------------------------------------------------

xi = dat$x
final_fit <- fit$fit[[fit$index$`min`[1]]][[fit$index$`min`[2]]]

fitcure <- final_fit$fitcure
predcure <- predict(fitcure, newx = xi,
                    s = min(fitcure$lambda), type = "response")

fitcox <- fi nal_fit$fitcox
predsurvexp <- predict(fitcox, newx = xi, s = min(fitcox$lambda), type = "response")


haz <- final_fit$haz
cumhaz <- final_fit$cumhaz

predsurv <- exp(-cumhaz*predsurvexp)

# Create Data frame to pass ----------------------------------------------

data = data.frame(dat$x) %>%
  mutate(times = dat$t, 
         event = dat$d,
         "truth_survival" = Surv(dat$t, dat$d))

# Predict Function (glmnet cure) ---------------------------------------------

# Trying to get predictions in similar format as tidymodels yardstick
# Eventually will integratethis with existing predict method function

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
    eval_timepoints = sorted_event_times
  }


  # Get Relative Cure and Cox Risk ------------------------------------------

  # Cure probability risk prediction (independent of timepoint)
  fitcure <- final_fit$fitcure
  predcure <- predict(fitcure, newx = data_x,
                      s = min(fitcure$lambda), type = "response")
  
  # Cox survival probability risk prediction (independent of timepoint)
  fitcox <- final_fit$fitcox
  predsurvexp <- predict(fitcox, 
                         newx = data_x,
                         s = min(fitcox$lambda),
                         type = "response")


  
  # Get Predicted probabilities at select timepoints -------------------
  
  # vector for conditional survival prob for all patients at sorted_survival_times[l]
  all_preds_all_tps <- list()
  
  # for each eval timepoint, get predicted survival values
  for (l in 1:length(eval_timepoints)) {
    all_preds <- data.frame("eval_timepoint" = NA, 
                            "preds" = NA, 
                            "pred_cure" = NA, 
                            "preds_surv" = NA, 
                            "weights" = NA)
    
    # vector for conditional survival prob for all patients at sorted_survival_times[l]
    predsurv <- rep(NA,length(survival_times)) 
    

    if (eval_timepoints[l] >= min(sorted_event_times)){
      
      ids <- which(
        sorted_event_times == 
          max(sorted_event_times[sorted_event_times <= eval_timepoints[l]])
        )
      
      # Conditional survival prob at given patient
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
      
      "preds" = preds, 
      "predcure" = predcure, 
      "predsurv" = predsurv, 
      "weights" = ipw) %>%
      mutate(id_num = 1:nrow(.))
    
    names(all_preds) <- c(
                          "preds",
                          "predcure", "predsurv", "weights",
                          "id_num")
    
    x <-  all_preds %>% nest(.preds = -c(id_num))
    
    all_preds_all_tps[[l]] <- all_preds
                            
  #  names(all_preds) = c("pred", "pred_cure","pred_surv")
  }
  # all_pred_df <- tibble(
  #   times = data$times,
  #   events = data$event,
  #   .pred = all_preds_all_tps
  # )
  return(all_pred_df)
}

predict_df <- predict_cure(final_fit = final_fit, data = data)


# Brier Calculation -------------------------------------------------------

brier_cure <- function(predict_df, times, event) {
  
  eval_timepoints = unique(predict_df$eval_timepoint)
  
  for (l in 1:length(eval_timepoints)) {
    brier[l] <- mean(I(di==0)*(1 - preds)^2/(ipw+0.001) + I(di==1)*(0 - preds)^2/(ipw+0.001))

  
  }
  
