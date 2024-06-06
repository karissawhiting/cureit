library(tidymodels)
library(censored)
tidymodels_prefer()


# Try Tidymodels Example --------------------------------------------------


# * Fit Lung Model  --------

data(cancer, package="survival")
data(lung_surv)

#lung <- lung %>% drop_na()

lung <- lung %>%
  mutate(status = case_when(
    status == 1 ~ 0, 
    status == 2 ~ 1
  ))


ph_spec <- 
  proportional_hazards() %>%
  set_engine("survival") %>% 
  set_mode("censored regression") 

ph_spec

set.seed(1)
ph_fit <- ph_spec %>% fit(Surv(time, status) ~ ., data = lung)
ph_fit

cox <- coxph(Surv(time, status) ~ ., data = lung)
cox
cumhaz_all<- survfit(cox)
cumhaz <- bind_cols(time = cumhaz_all$time, surv = cumhaz_all$surv)

#cumhaz$surv[which.min(abs(cumhaz$time - closest_event_time))]

#cumhaz$surv[cumhaz$time[]

# * Predict Lung Model  -----------

# Make predictions
p <- predict(
  ph_fit, 
  lung, 
  type = "survival",
  eval_time = lung$time
) %>%
  mutate(surv_truth =  Surv(lung$time, lung$status))

# add weights
p_weights <- .censoring_weights_graf(ph_fit, predictions = p) %>%
  mutate(pat_id = 1:nrow(.))

p_weights_long <- p_weights %>%
  unnest(cols = .pred)

# calc brier
brier_scores_lung <-
  p_weights %>% 
  brier_survival(truth = surv_truth, .pred)

brier_scores_lung

# Try our code on their data ----------------------------------------------

# 
# predict_df_long <- predict_df %>%
#   unnest(everything())

predict_df_nest_time <- p_weights_long %>%
  nest(data = -.eval_time)

predict_df_timepoint = predict_df_nest_time[10,] %>%
  unnest(everything())

brier_for_each_eval <- function(predict_df_timepoint) {
  
  
  x = mean(I(predict_df_timepoint$times > predict_df_timepoint$.eval_time )*(1 - predict_df_timepoint$.pred_survival)^2/(predict_df_timepoint$.weight_censored+0.001) + 
             I(predict_df_timepoint$event == 1 & predict_df_timepoint$times <= predict_df_timepoint$.eval_time  )*
             (0 - predict_df_timepoint$.pred_survival)^2/(predict_df_timepoint$.weight_censored+0.001))
  x
}




lung2 <- lung %>%
  rename("times" = time, 
         "event" = status)
cox <- coxph(Surv(times, event) ~ ., data = lung2)
cox
cumhaz_all<- survfit(cox)
cumhaz <- bind_cols(time = cumhaz_all$time, surv = cumhaz_all$surv)


lung2 <- lung2 %>%
  mutate(truth_survival = Surv(times, event))

predict_cure(final_fit = cox, new_data = lung2, eval_timepoints = lung2$times)




# Misc --------------------------------------------------------------------


eval_time_test <- 300

new_data = lung
survival_times <- lung$time
sorted_survival_times <- sort(survival_times)

survival_status <- lung$status
sorted_event_times = sort(survival_times[survival_status==1])

censoring_survfit = survfit(Surv(survival_times, 1 - survival_status) ~ 1)
censoring_weights = summary(censoring_survfit, times = survival_times)
censoring_weights = censoring_weights$surv

all_weights = bind_cols(censoring_weights_time, censoring_weights)
names(all_weights) <- c("weights_time", "weights_surv")

eval_timepoint = 300

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
  censoring_weight_at_eval_time <- censoring_weights[timepoint_index]
  
  ipw <- 
    as.numeric(survival_times > eval_timepoint)*censoring_weight_at_eval_time + 
    as.numeric(survival_times <= eval_timepoint)*censoring_weights
  
  ipw = 1/ipw
  
  # if no events occurred before eval time 
} else if (eval_timepoint < min(sorted_event_times)) {
  
  predsurv <- rep(1,length(survival_times))
  ipw <- 1
  
}

cbind(survival_times, ipw) %>% View()

# Compare ------
lung_surv_long <- lung_surv %>% unnest(everything())

lung_surv_long <- lung_surv_long %>% filter(.eval_time == eval_time_test)

#lung %>% select(time, status) %>% distinct() %>% nrow()


# Here- New ---------------------------------------------------------------------
# ti & di is censored patient - If censored time is after eval time () ti > tbrier[l]
# if event, and event time earlier than evaluation time - (0 - preds)^2/(ipw+0.001))

# 1) If the observed time is a censoring time and that is before the evaluation time, the data point should make no contribution to the performance metric (their "category 3"). These values have a missing value for their probability estimate (and also for their weight column).
# ^ not included, ensure they are NAs in predict df 

# 2)  I(ti <= tbrier[l] & di == 1)*(0 - preds)^2/(ipw+0.001)) 
# Old version ----
# brier[l] <- mean(I(ti > tbrier[l] & di==0)*
#                    (1 - preds)^2/(ipw+0.001) + 
#                    
#                    I(ti <= tbrier[l] & di == 1)*(0 - preds)^2/(ipw+0.001))

# Updated Version ---


brier[l] <- mean(I(ti > tbrier[l])*
                   (1 - preds)^2/(ipw+0.001) + 
                   
                   I(ti <= tbrier[l] & di == 1)*(0 - preds)^2/(ipw+0.001))

