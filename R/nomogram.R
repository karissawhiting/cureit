# NOMOGRAM SCRIPT -------------------------------------------------------------
# -----------------------------------------------------------------------------

library(tidyverse)

# fit an smcure model but with terms and data attached for multi-level factors
cure_nomogram <- function(fit, prediction_time = 3) {
  

  # extract from model fit---
  var_metadata <- fit$variable_metadata
  fit_obj <- fit$cure_fit_obj
  mm <- fit$model_matrix
  betas <- fit_obj$beta
  data <- fit$data
  
  # change to be only if backticks
  names(betas) <- str_remove(names(betas), fixed("X[, -1]")) %>%
    str_remove_all(., fixed("`"))
  
  
  # separate continuous and categorical ---
  cont <- var_metadata %>% 
    filter(factor == "no") %>%
    pull(original_vars)
  
   cont2 <- var_metadata %>% 
    filter(factor == "no") %>%
    pull(new_var)
   
  cat2 <- var_metadata %>% 
    filter(factor == "yes") %>%
    pull(new_var)
  
   beta_cont <- betas[names(betas) %in% cont2]
   beta_cat <- betas[names(betas) %in% cat2]
   
  cont_df <- data %>%
    dplyr::select(all_of(cont))

  
  # make "pretty" continuous ---
   cont_pretty <- map(cont_df, ~pretty(.x, 10))
   cont_pretty_beta <-   map2(cont_pretty, beta_cont, 
                              ~.x*.y)
  
     
  # make pretty categorical ---
  # cat <- var_metadata %>% 
  #   filter(factor == "yes") %>%
  #   pull(original_vars) %>%
  #   unique()
  # 
  # cat_df <- data %>%
  #   select(all_of(cat))
  
   cat_df <- mm %>%
     tibble::as_tibble() %>%
     select(all_of(cat2))
   
  cat_pretty <- map(cat_df,
                    ~as.factor(rep(levels(as.factor(.x)), length.out = 10)))
  
#  cat_pretty_bin <- map(cat_df, ~rep(as.numeric(as.factor(levels(as.factor(.x)))), length.out = 10))
    
  
    cat_pretty_beta <-   map2(cat_df, beta_cat, 
                              ~unique(.x*.y))
   
    all_pretty_beta <- c(cont_pretty_beta, cat_pretty_beta)
    
  # calculate x*Beta ---
  x_beta <- t(t(mm) * fit_obj$beta)
  
  
  
  # get scaling constant ------------------------------------------------------
   # ranges <- apply(x_beta, 2, function(x) {
   #   m1 <- max(x) 
   #   m2 <- min(x)
   #   m1 - m2})
   # 

    ranges <- map_dbl(all_pretty_beta, function(x) {
     m1 <- max(x) 
     m2 <- min(x)
     m1 - m2})
   
   # for use later 
  mins <- map_dbl(all_pretty_beta, min)
 # max <- map_dbl(all_pretty_beta, max)
  stats <- cbind.data.frame(mins, ranges) %>%
    tibble::rownames_to_column(var = "new_var") %>%
    left_join(select(var_metadata, original_vars, new_var)) %>%
    group_by(original_vars) %>%
    mutate(variable_min = min(mins))
  
  max_range <- max(ranges)
   sc <- 100/max_range

  
  # COX NOMOGRAM --------------------------------------------------


  # Coefficients ---
  
  
#   coefs <- apply(x_beta, 2, function(x) (x - min(x))*sc) %>%
 #    tibble::as_tibble()
  
  x_beta_tb <- x_beta %>%
    tibble::as_tibble()
   
  names(all_pretty_beta) <- str_remove(names(all_pretty_beta), fixed("X[, -1]")) %>%
    str_remove_all(., fixed("`"))
  
#  all_pretty_beta_0 <- map(all_pretty_beta, ~.x[.x > 0])
  
  coefs <- map2(all_pretty_beta, 
                   stats$variable_min,
                   function(x, y) (x - y)*sc)
  

  cat_final_vars  <- tibble(coefs) %>% bind_cols(names(all_pretty_beta)) %>%
    tidyr::unnest(cols = c(coefs)) %>%
    rename("new_var" = `...2`) %>%
    left_join(select(var_metadata, original_vars, new_var, factor)) %>%
    filter(factor == "yes")
 
  
  # x <- all_pretty_beta[[1]]
  # x - 1
  # 
  #  df_coefs <- as.data.frame(coefs) %>%
  #    mutate(across(everything(), ~round(.x, 1)))
  #  
  #  # categorical
  # cat <- var_metadata %>% 
  #   filter(factor == "yes") %>%
  #   pull(new_var)
  #   
  #  df_coefs_cat <- df_coefs %>% 
  #    select(all_of(cat))
  #  
  # final_cat_vars <- df_coefs_cat %>%
  #   pivot_longer(everything()) %>%
  #   left_join(., var_metadata, by = c("name" = "new_var")) %>%
  #   mutate(display_level = case_when(value == 0 ~ ref_level, 
  #                                    TRUE ~ name)) %>%
  #   select(display_level, name, value, original_vars) %>%
  #   transmute(x = value, 
  #          label = display_level, 
  #          y = original_vars) %>%
  #   mutate(label = str_remove(label, y))
   
  final_cat_vars  <-  cat_final_vars %>%
    left_join(., var_metadata) %>%
    mutate(display_level = case_when(coefs == 0 ~ ref_level, 
                                     TRUE ~ new_var)) %>%
    select(display_level, new_var, coefs, original_vars) %>%
    transmute(x = coefs, 
           label = display_level, 
           y = original_vars) %>%
    mutate(label = str_remove(label, y))

  
  final_cont_vars <- pmap_df(list(beta_cont, cont_pretty, names(cont_pretty)),
                            function(x, y, z) {
    cbind.data.frame(x = round(x*y, 2), label = as.character(y)) %>%
      mutate(y = z) })
    
    
  # Points  ----------
  # Linear Predictor ---
  
  lp_data <- x_beta 
  lp_data <- cbind.data.frame(lp_data, lp_raw = rowSums(lp_data)) %>%
    mutate(lp = scale(lp_raw, scale = FALSE))
  
  Intercept_raw <- -mean(lp_data$lp_raw)
  Intercept <- Intercept_raw  + sum(mins)
  
  lp <- pretty(range(lp_data$lp), n = 10)
  scaled_lp <- (lp-Intercept) *sc
  
  
  points = pretty(0:100, 10)
  df_points <- as.data.frame(points) %>%
    transmute(x = points,
      label = as.character(points), 
      y = "Points")
  
  
  upper_range <-  max(sc*max(lp-Intercept), 100)
  total_points <- pretty(c(0, upper_range), n=10)
  upper_range_pretty <- max(total_points)
  
  df_lp <- as.data.frame(lp) %>%
    transmute(x = scaled_lp, 
      label = as.character(lp), 
      y = "Linear Predictor") %>%
    mutate(x = (x/upper_range_pretty)*100)
  
  
  df_tp <- as.data.frame(total_points) %>%
    mutate(label = as.character(round(total_points, 0)), 
      x = scales::rescale(total_points, to = c(0, 100)), 
      y = "Total Points") %>%
    select(-total_points)
  
    
  # Prediction -------
  
  baseline <- cbind.data.frame(s = fit_obj$s, time = fit_obj$Time) %>%
    mutate(diff = abs(time- prediction_time)) %>%
    filter(diff == min(abs(diff)))
  
  s0 <- baseline$s[1]
  
  df_pred1 <- as.data.frame(scaled_lp) %>%
    transmute(x = (scaled_lp/upper_range_pretty)*100, 
      label = as.character(round(s0^exp((scaled_lp/sc) + Intercept + 
          mean(lp_data$lp_raw)), 1)), 
      y = paste0(prediction_time,  " Surv Prob"))

     
  # Plot -----------------
  
  all <- bind_rows(final_cat_vars, final_cont_vars, df_points, df_pred1)
  
  all <- all %>%
    mutate(y = fct_relevel(y, "Points")) %>%
    mutate(y = fct_rev(y))
  
  p1 <- all %>%
    ggplot(aes(x = x, y = y)) + geom_line() +
    geom_point() + 
    geom_text(aes(label = label), vjust = 1.5)  + ylab(" ") + xlab(" ") +
    ggtitle("Estimated Survival for Uncured") +
    theme_minimal() +
    theme(panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.line = element_blank())
  
  p1
  
}

