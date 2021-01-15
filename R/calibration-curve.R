

### Calibration Curve -----------------------------------------------------------

cure_calibration <- function(fit, prediction_time = 3) {
  fit_obj <- fit$smcure_model_object
  mm <- fit$model_matrix
  data <- fit$data
  var_metadata <- fit$variable_metadata
  formula <- fit$formula
  mf <- formula
  mf <- eval(mf, parent.frame())
  
  p <- predictsmcure(fit_obj,
    newX = mm,
    newZ = mm
  )


  preds <- as_tibble(p$prediction) %>%
    #  bind_cols(Time = fit$Time) %>%
    mutate(time_dif = abs(Time - prediction_time)) %>%
    filter(time_dif == min(time_dif)) %>%
    dplyr::select(-Time, -time_dif) %>%
    t() %>%
    as_tibble() %>%
    mutate(index = 1:nrow(.)) %>%
    mutate(decile = ntile(V1, 10))


  data <- data %>%
    mutate(index = 1:nrow(.))

  mean_pred <- preds %>%
    group_by(decile) %>%
    summarise(mean_pred = mean(V1))

  one <- preds %>%
    filter(decile == 2)

   
  survfit_form <- as.formula(paste0("Surv(", 
                         all.vars(mf)[1], " , ",
                          all.vars(mf)[2], 
                         ") ~ 1"))
  
  f <- function(one) {
    x <- data %>% filter(index %in% one$index)
    t <- survfit(survfit_form,
      data = x
    )
    broom::tidy(t) %>%
      mutate(time_dif = abs(time - prediction_time)) %>%
      filter(time_dif == min(time_dif)) %>%
      pull(estimate)
  }


  preds2 <- preds %>%
    group_by(decile) %>%
    nest() %>%
    mutate(obs_prob = map_dbl(data, ~ f(.x)))


  all <- preds2 %>%
    left_join(mean_pred) %>%
    dplyr::select(-data)

  all %>%
    ggplot(aes(x = mean_pred, y = obs_prob)) +
    geom_point() +
    geom_line() +
    theme_minimal() +
    geom_abline(slope = 1, linetype = 2) 
}

# Example 

# fit <- fit_cure(
#   formula = Surv(days, status) ~ ulc + sex + thick_cat,
#   data = mel
# )

# cure_calibration(fit, prediction_time =300)
