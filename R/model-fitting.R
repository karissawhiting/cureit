# Functions to Run Cure Models (Once or Multiple Times) ------------------------
# ------------------------------------------------------------------------------


# Single Model Run -------------------------------------------------------------
# Wrapper function for smcure to accept multi-level variables and save necessary 
# Data for Nomogram

fit_cure <-  function(formula, data) { 
  
  # pull original variables from formula
   mf <- formula
   mf <- eval(mf, parent.frame())
   
   original_var_names <- bind_cols(labels(terms(mf)),
                               1:length(labels(terms(mf)))) %>%
     set_names("original_vars", "index")
   
   # only keep variables we need 
    data_sub <- data %>%
      select(all.vars(mf)[1],
            all.vars(mf)[2], 
            labels(terms(mf))) %>%
    mutate(across(is.character, ~as.factor(.x)))
   
   # create model matrix for dummy variables
   model_mat <- model.matrix(mf, data = data_sub)
   model_mat_sub <- model_mat[,-1]
   
   # edited data set for necessary variables in smcure
  data_for_smcure <- data_sub %>%
      select(all.vars(mf)[1],
            all.vars(mf)[2])%>%
     bind_cols(as.data.frame(model_mat_sub))
   
   # create new formula from dummy variables for cure model
   vars <- colnames(model_mat_sub)
   w <- grep(pattern = " ", x = vars, fixed = TRUE)
   vars[w] <- sprintf("`%s`", vars[w])
   
   new_formula <- paste0("Surv(", all.vars(mf)[1], ", ", all.vars( mf)[2], ")",
                         " ~ ", paste(vars, collapse=" + ") )
   surv_formula <- as.formula(new_formula)
   
   # NOTE: eventually need to allow different cure variables 
   cure_formula <- surv_formula[-2]
   
   # results -----
   # [1] run model (add other variables here)
   model <- smcure(surv_formula,
          cureform = cure_formula,
       data = data_for_smcure,
    model="ph", Var = TRUE)

   # [1] create variable metadat
   a <- attributes(model_mat)
   levels <- map(data_sub, ~levels(.x)) %>%
     discard(is.null)
   
   var_metadata <- bind_cols(a$dimnames[[2]],
                             a$assign) %>% 
     set_names("new_var", "index") %>%
     left_join(original_var_names) %>%
     filter(new_var != "(Intercept)") 
   
   var_metadata <- var_metadata %>%
     mutate(levels = map(original_vars, ~levels(data_sub[[.x]]))) %>%
     mutate(other_levels = str_remove(new_var, original_vars)) %>%
     group_by(original_vars) %>% 
     group_by(original_vars) %>%
     mutate(other_levels = list(other_levels)) %>%
     mutate(ref_level = map2(levels, other_levels, ~unique(setdiff(.x, .y)))) %>%
     group_by(index) %>%
     mutate(n = n()) %>%
     mutate(factor = case_when(original_vars %in% names(a$contrasts) ~ "yes", 
                               TRUE ~ "no")) %>%
     ungroup() %>%
     unnest(., ref_level, keep_empty = TRUE)
   
  
   return(list("smcure_model_object" = model,
               "variable_metadata" = var_metadata,
               "model_matrix" = model_mat_sub, 
               "data" = data_sub, 
               "formula" = formula))
}


# Multiple Model Runs -----------------------------------------------------------

# runs model a number of times and returns mode of p-values 
multiple_mod_runs <- function(formula,
                              nboot = 200,
                              eps = 0.0001, 
                              num_repeats = 10, 
                              data = mel) {

  source(here::here("R", "clean-model-results.R"))
  
  # pull original variables from formula
  mf <- formula
  mf <- eval(mf, parent.frame())

  # original_var_names <- bind_cols(labels(terms(mf)),
  #                             1:length(labels(terms(mf)))) %>%
  #   set_names("original_vars", "index")
  #
  # only keep variables we need
  data_sub <- data %>%
    select(
      all.vars(mf)[1],
      all.vars(mf)[2],
      labels(terms(mf))
    ) %>%
    mutate(across(is.character, ~ as.factor(.x)))

  # create model matrix for dummy variables
  model_mat <- model.matrix(mf, data = data_sub)
  model_mat_sub <- model_mat[, -1]

  # edited data set for necessary variables in smcure
  data_for_smcure <- data_sub %>%
    select(
      all.vars(mf)[1],
      all.vars(mf)[2]
    ) %>%
    bind_cols(as.data.frame(model_mat_sub))

  # create new formula from dummy variables for cure model
  vars <- colnames(model_mat_sub)
  w <- grep(pattern = " ", x = vars, fixed = TRUE)
  vars[w] <- sprintf("`%s`", vars[w])

  new_formula <- paste0(
    "Surv(", all.vars(mf)[1], ", ", all.vars(mf)[2], ")",
    " ~ ", paste(vars, collapse = " + ")
  )
  surv_formula <- as.formula(new_formula)

  # NOTE: eventually need to allow different cure variables
  cure_formula <- surv_formula[-2]

  # model runs
  repeated_model_run <- function() {
    tryCatch(
      {
        smcure(surv_formula,
          cureform = cure_formula,
          data = data_for_smcure,
          model = "ph", Var = T,
          nboot = nboot,
          eps = eps
        )
      },
      error = function(e) {
        NA
      }
    )
  }
  
  res <- map(1:num_repeats, ~repeated_model_run())
  
  res <- res[!is.na(res)]
  
  message(paste0(length(res), " of ", num_repeats, " successful model runs"))
  
  res_cure <- map_dfr(res, ~mmcure(.x)) %>%
    filter(name != "(Intercept)") 
  
  var_cure_stab <- res_cure %>%
    ggplot(aes(x = value, fill = name)) + geom_density(alpha = .8)
  
  p_cure_stab <- res_cure %>%
    ggplot(aes(x = pvalue, fill = name)) + geom_density(alpha = .8)
  
   den_fun <- function(dd) {
     final_p <- dd$x[which.max(dd$y)]
     final_p
   }
   
  
   
 cure_site <- res_cure %>%
   group_by(name) %>%
   nest() %>%
   mutate(density = map(data, ~density(.x$pvalue))) %>%
   mutate(final_p = map_dbl(density, ~den_fun(.x))) %>%
   mutate(est = map_dbl(data,
                        ~.x$value[which.min(abs(round(.x$pvalue, 2) - round(final_p, 2)))] %>%
    unique())) %>%
   transmute(variable = str_remove(name, fixed("Z[, -1]")), 
             value_cure = est, 
             p_value_cure = final_p) %>%
   mutate(variable = str_remove_all(variable, fixed("`"))) %>%
      ungroup() %>%
   dplyr::select(-name)
          
 # Surv-----------
  
  res_surv <- map_dfr(res, ~mmsurv(.x))
 
   var_surv_stab <- res_surv %>%
    ggplot(aes(x = value, fill = name)) + geom_density(alpha = .8)
  
  p_surv_stab <- res_surv %>%
    ggplot(aes(x = pvalue, fill = name)) + geom_density(alpha = .8)
    
 surv_site <- res_surv %>%
   group_by(name) %>%
   nest() %>%
   mutate(density = map(data, ~density(.x$pvalue))) %>%
   mutate(final_p = map_dbl(density, ~den_fun(.x))) %>%
   mutate(est = map_dbl(data,
                        ~.x$value[.x$pvalue == round(final_p, 2)] %>%
    unique())) %>%
   ungroup() %>%
   transmute(variable = str_remove(name, fixed("X[, -1]")), 
             value_surv = est, 
             p_value_surv = final_p) %>%
      mutate(variable = str_remove_all(variable, fixed("`")))
 
  all_models <- left_join(cure_site, surv_site)


  all_models <- all_models %>%
    mutate(variable = 
             tools::toTitleCase(str_replace_all(variable, "_", " "))) %>%
    mutate_if(is.numeric, ~round(.x, 2))
  
  return(list("model_results" = all_models, 
              "var_surv_stab" = var_surv_stab, 
              "p_surv_stab" = p_surv_stab, 
              "var_cure_stab" = var_cure_stab, 
              "p_cure_stab" = p_cure_stab))

}


  