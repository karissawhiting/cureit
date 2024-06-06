#' Predict Cure and Survival Probabilities
#'
#' This function predicts cure and survival probabilities for data based on a final fitted model.
#' It calculates predicted survival probabilities at specified evaluation time points, accounting for censoring.
#' Output is consistent with tidymodels predict guidelines and can be passed into yardstick functions like `brier()`
#'
#' @param final_fit A list containing the fitted cure and survival models. It should have elements `fitcure`, `fitcox`, `haz`, and `cumhaz`.
#' @param new_data A data frame containing the new observations to predict on. It should include columns `times` (observed survival times), `event` (event indicator, 1 for event, 0 for censoring), and other covariates used for prediction.
#' @param eval_timepoints Optional. A vector of time points at which to evaluate survival probabilities. If NULL, the function uses the observed event times from `new_data`.
#' 
#' @return A data frame containing the predicted survival probabilities, cure rates, and other related metrics for each observation at the specified evaluation time points. The returned data frame includes the following columns:
#' \itemize{
#'   \item \code{id}: Identifier for each observation.
#'   \item \code{.eval_time}: Evaluation time point.
#'   \item \code{.pred_survival}: Predicted survival probability at the evaluation time point.
#'   \item \code{.pred_curerate}: Predicted cure rate for each observation.
#'   \item \code{.pred_survportion}: Component of the survival prediction that excludes the cure component.
#'   \item \code{.weight_censored_raw}: Raw censoring weight at each evaluation time point.
#'   \item \code{.weight_censored}: Inverse probability weight for censoring.
#'   \item \code{event}: Event indicator from the new data.
#'   \item \code{times}: Observed survival times from the new data.
#' }
#'
#' @examples
#' \dontrun{
#'   # Assuming `final_fit` is a fitted model list and `new_data` is a data frame with appropriate structure:
#'   predict_df <- predict_cure(final_fit = final_fit, new_data = new_data)
#'   print(predict_df)
#' }
#' 
#' @importFrom dplyr select mutate bind_cols
#' @importFrom tidyr nest unnest
#' @importFrom tibble tibble
#' 
#' @export
predict_cure <- function(final_fit,
                         new_data = NULL,
                         eval_timepoints = NULL) {
  
  # Values from Observed Data & Model Fit ---------------------------------
  
  survival_times = new_data$times
  events = new_data$event
  sorted_survival_times <- sort(survival_times)
  sorted_event_times = sort(new_data$times[new_data$event==1])
  
  # Extract cumulative hazard from final fit (for event times)
  haz <- final_fit$haz
  cumhaz <- final_fit$cumhaz
  
  # Set up Data 
  data_x <- new_data %>%
    select(-c(times, event, truth_survival)) %>%
    as.matrix()
  
  # Evaluation Time Points --------------------------------------------------
  # Default evaluation time points are observed event times if not supplied by user
  
  all_eval_timepoints <- eval_timepoints %||% sorted_survival_times
  
  # Censoring Weights --------------------------------------------------------
  
  censoring_survfit = survfit(Surv(survival_times, 1 - events) ~ 1)
  
  # CHECK- or censoring weights at all all respective survival times?
  # censoring_weights = summary(censoring_survfit, times = all_eval_timepoints)
  censoring_weights = summary(censoring_survfit, times = survival_times)$surv

    
  # Relative Risk for Cure and Cox Risk ------------------------------------
  # (independent of time point) 
  
  #  * Cure probability risk prediction ------
  fitcure <- final_fit$fitcure
  predcure <- predict(fitcure, newx = data_x,
                      s = min(fitcure$lambda), type = "response")
  
  # * Cox survival probability risk prediction ------
  fitcox <- final_fit$fitcox
  predsurvexp <- predict(fitcox, 
                         newx = data_x,
                         s = min(fitcox$lambda),
                         type = "response")
  
  
  
  # Predicted Probabilities at Evaluated Time points -------------------------
  # (dependent on time point- need IPCW)
  
  all_preds_all_tps <- list()
  
  # Evaluated in all patients, for each eval time point (iteration on time point)
  for (timepoint_index in 1:length(all_eval_timepoints)) {
    
    # timepoint we are evaluating
    eval_timepoint <- all_eval_timepoints[timepoint_index] 
    

    # * Inverse Probability Censoring Weights ---------
    
    # If no time points in new data less than lowest event time in model fit data (most cases):
    if (eval_timepoint >= min(sorted_event_times)) {
      
      # get maximum event time that is less than time point we are evaluating 
      closest_event_time <- max(
        sorted_event_times[sorted_event_times <= eval_timepoint])
      
      index_closest_event_time <- which(sorted_event_times == closest_event_time)
      
      # Get predicted survival by using cumulative hazard at closest previous event time
      predsurv <- exp(-cumhaz[index_closest_event_time]*predsurvexp)
      
      # CHECK - Probability of the censoring at evaluation time point (or should it be at previous index closest event time??)
     #  censoring_weight_at_eval_time <- censoring_weights[timepoint_index]
      censoring_weight_at_eval_time <- censoring_weights[index_closest_event_time]
      
      ipw <- 
        # if observed time is greater than event time, use censoring weight at evaluation time (or time closest)
        (as.numeric(survival_times > eval_timepoint)*censoring_weight_at_eval_time) + 
        
        # if observed time is less than or equal to eval time, use censoring weight at it's own observed time
        as.numeric(survival_times <= eval_timepoint)*censoring_weights
      
      ipw = 1/(ipw+ 0.001)
      
      # if no events occurred before eval time 
    } else if (eval_timepoint < min(sorted_event_times)) {
      
      predsurv <- rep(1,length(survival_times))
      ipw <- 1
      
    }
    
    preds <- 1 - predcure + predcure*predsurv
    
    all_preds <- tibble(
      
      ".eval_time" = eval_timepoint,
      ".pred_survival" = preds, 
      ".pred_curerate" = predcure, 
      ".pred_survportion" = predsurv, 
      ".weight_censored_raw" = censoring_weights[timepoint_index],
      ".weight_censored" = ipw) %>%
      mutate(id = 1:nrow(.))
    
    all_preds_nest <-  all_preds %>% nest(.preds = -c(id))
    
    all_preds_all_tps[[timepoint_index]] <- all_preds
    
  }
  
  all_pred_df <- all_preds_all_tps %>% 
    do.call("rbind", .)
  
  all_pred_df <- all_pred_df %>%
    nest(.pred = -c(id)) %>% 
    bind_cols(., "event" = new_data$event, 
              "times" = new_data$times)
  
  return(all_pred_df)
  
}


# Try Function ------------------------------------------------------------
# 
# predict_df <- predict_cure(final_fit = final_fit, new_data = data)
# predict_df_long <- predict_df %>%
#   unnest(everything())


