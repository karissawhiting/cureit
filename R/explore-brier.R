library(tidymodels)
library(censored)
tidymodels_prefer()


# Try Tidymodels Example --------------------------------------------------


# Fit Lung Model  ---------------------------------------------------------

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
p_weights <- .censoring_weights_graf(ph_fit, predictions = p)

# calc brier
brier_scores_lung <-
  p_weights %>% 
  brier_survival(truth = surv_truth, .pred)


# Try our code on their data ----------------------------------------------


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
