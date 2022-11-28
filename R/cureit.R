#' Cure Model Regression
#'
#' @param object input object
#' @param surv_formula formula with `Surv()` on LHS and covariates on RHS.
#' @param cure_formula formula with covariates for cure fraction on RHS
#' @param data data frame
#' @param conf.level confidence level. Default is 0.95.
#' @param nboot number of bootstrap samples used for inference.
#' @param eps convergence criterion for the EM algorithm.
#' @param ... passed to methods
#'
#' @return cureit object
#' @family cureit() functions
#' @name cureit
#' @examples
#' cureit(surv_formula = Surv(ttdeath, death) ~ age + grade, 
#' cure_formula = ~ age + grade,  data = trial)
NULL

# Formula method
#' @rdname cureit
#' @export
cureit.formula <- function(surv_formula, cure_formula, data, conf.level = 0.95, nboot = 100, eps = 1e-7,...) {
  
  # checking inputs  -------------------------
  checking(surv_formula = surv_formula, cure_formula = cure_formula, data = data)
  
  # process model variables ----------------------------------------------------
  processed <- cureit_mold(surv_formula, cure_formula, data)
  
  # building model -------------------------------------------------------------
  cureit_bridge(processed, surv_formula, cure_formula, data, conf.level = conf.level, nboot = nboot, eps = eps)
}

cureit_mold <- function(surv_formula, cure_formula, data, surv_blueprint = NULL, cure_blueprint = NULL) {
  
  if (is.null(surv_blueprint)) surv_blueprint = hardhat::default_formula_blueprint(intercept = TRUE)
  
  surv_processed <-
    hardhat::mold(
      surv_formula, data,
      blueprint = surv_blueprint
    )
  # remove intercept
  surv_processed$predictors <- surv_processed$predictors[, -1]
  surv_processed$predictors <- surv_processed$predictors %>% janitor::clean_names()
  surv_processed$outcomes <- surv_processed$outcomes %>% janitor::clean_names()
  
  
  if (is.null(cure_blueprint)) cure_blueprint = hardhat::default_formula_blueprint(intercept = TRUE)
  
  cure_processed <-
    hardhat::mold(
      cure_formula, data,
      blueprint = cure_blueprint
    )
  # remove intercept
  cure_processed$predictors <- cure_processed$predictors[, -1]
  cure_processed$predictors <- cure_processed$predictors %>% janitor::clean_names()
  
  list(surv_processed=surv_processed,cure_processed=cure_processed)
}

checking <- function(surv_formula, cure_formula, data, keep_all = FALSE) {
  # evaluating LHS of formula --------------------------------------------------
  surv_formula_lhs <-
    tryCatch(
      {
        rlang::f_lhs(surv_formula) %>%
          rlang::eval_tidy(data = data)
      },
      error = function(e) {
        cli::cli_alert_danger("There was an error evaluating the LHS of the formula.")
        stop(e, call. = FALSE)
      }
    )
  
  ### !!! Need to check cure model formula as well: LHS should be empty
  
  # checking type of survival model formula LHS -------------------------------------------------------
  if (!inherits(surv_formula_lhs, "Surv") ||
      !identical(attr(surv_formula_lhs, "type"), "right")) {
    paste(
      "The LHS of the formula must be of class 'Surv' and type 'right'.",
      "Please review expected syntax in the help file.",
      "The status variable must be a factor, where the first level indicates",
      "the observation was censored, and subsequent levels are the",
      "competing events. Cannot use `Surv(time2=)` argument."
    ) %>%
      stop(call. = FALSE)
  }
  
}

new_cureit <- function(surv_coefs, surv_coef_names, cure_coefs, cure_coef_names, 
                       surv_formula_input, cure_formula_input, surv_formula_smcure,
                       cure_formula_smcure, tidy, smcure, data,
                       surv_blueprint, cure_blueprint, conf.level, nboot, eps) {
  
  # function to create an object
  
  if (!is.numeric(surv_coefs)) {
    stop("`surv_coefs` should be a numeric vector.", call. = FALSE)
  }
  
  if (!is.character(surv_coef_names)) {
    stop("`surv_coef_names` should be a character vector.", call. = FALSE)
  }
  
  if (length(surv_coefs) != length(surv_coef_names)) {
    stop("`surv_coefs` and `surv_coef_names` must have the same length.")
  }
  
  if (!is.numeric(cure_coefs)) {
    stop("`cure_coefs` should be a numeric vector.", call. = FALSE)
  }
  
  if (!is.character(cure_coef_names)) {
    stop("`cure_coef_names` should be a character vector.", call. = FALSE)
  }
  
  if (length(cure_coefs) != length(cure_coef_names)) {
    stop("`cure_coefs` and `cure_coef_names` must have the same length.")
  }
  
  hardhat::new_model(
    surv_coefs = surv_coefs %>% stats::setNames(surv_coef_names),
    cure_coefs = cure_coefs %>% stats::setNames(cure_coef_names),
    surv_formula = surv_formula_input,
    cure_formula = cure_formula_input,
    data = data,
    conf.level = conf.level,
    nboot = nboot,
    eps = eps,
    surv_xlevels =
      stats::model.frame(surv_formula_input, data = data)[, -1, drop = FALSE] %>%
      purrr::map(
        function(.x) {
          if (inherits(.x, "factor")) {
            return(levels(.x))
          }
          if (inherits(.x, "character")) {
            return(unique(.x) %>% sort())
          }
          return(NULL)
        }
      ) %>%
      purrr::compact(),
    cure_xlevels =
      
      stats::model.frame(cure_formula_input, data = data)[, , drop = FALSE] %>%
      
      purrr::map(
        function(.x) {
          if (inherits(.x, "factor")) {
            return(levels(.x))
          }
          if (inherits(.x, "character")) {
            return(unique(.x) %>% sort())
          }
          return(NULL)
        }
      ) %>%
      purrr::compact(),
    tidy = tidy,
    smcure = smcure,
    surv_blueprint = surv_blueprint,
    cure_blueprint = cure_blueprint,
    class = "cureit"
  )
}

cureit_impl <- function(surv_formula, cure_formula, newdata, conf.level = conf.level, nboot = nboot, eps=eps) {
  
  # function to run cureit and summarize with tidy (implementation)
  quiet_smcure <- purrr::quietly(smcure::smcure)
  
  if (nboot == 0){
    
    cureit_fit <-
      quiet_smcure(formula=surv_formula,
                   cureform=cure_formula,
                   data=newdata,
                   model="ph",
                   eps=eps,
                   Var = FALSE
      )
    
    cureit_fit <- cureit_fit$result
    cureit_fit$b_sd <- cureit_fit$b_zvalue <- cureit_fit$b_pvalue <- NA
    cureit_fit$beta_sd <- cureit_fit$beta_zvalue <- cureit_fit$beta_pvalue <- NA
    
    # broom method can be constructed later 
    tidy <- broom::tidy(cureit_fit, conf.int = FALSE)
    
    s.coefs <- tidy$df_surv$estimate
    s.coef_names <- tidy$df_surv$term
    c.coefs <- tidy$df_cure$estimate
    c.coef_names <- tidy$df_cure$term
    
    list(
      surv_coefs = s.coefs,
      surv_coef_names = s.coef_names,
      cure_coefs = c.coefs,
      cure_coef_names = c.coef_names,
      tidy = tidy,
      smcure = cureit_fit
    )
    
    
  }else if (nboot > 0){
    cureit_fit <-
      quiet_smcure(formula=surv_formula,
                   cureform=cure_formula,
                   data=newdata,
                   model="ph",
                   nboot=nboot,
                   eps=eps
      )
    
    cureit_fit <- cureit_fit$result
    
    # broom method can be constructed later 
    tidy <- broom::tidy(cureit_fit, conf.int = TRUE, conf.level = conf.level)
    
    s.coefs <- tidy$df_surv$estimate
    s.coef_names <- tidy$df_surv$term
    c.coefs <- tidy$df_cure$estimate
    c.coef_names <- tidy$df_cure$term
    
    list(
      surv_coefs = s.coefs,
      surv_coef_names = s.coef_names,
      cure_coefs = c.coefs,
      cure_coef_names = c.coef_names,
      tidy = tidy,
      smcure = cureit_fit
    )
    
  }
  
  
  
  
}

cureit_bridge <- function(processed, surv_formula_old, cure_formula_old, data, conf.level, nboot, eps) {
  
  # function to connect object and implementation
  
  # validate_outcomes_are_univariate(processed$outcomes)
  
  s.predictors <- as.matrix(processed$surv_processed$predictors)
  s.outcomes <- as.matrix(processed$surv_processed$outcomes[, 1, drop = TRUE])
  surv_formula <- paste0("Surv(", paste(colnames(s.outcomes),collapse=","), ")",
                         " ~ ", paste(colnames(s.predictors), collapse=" + ") )
  surv_formula <- as.formula(surv_formula)
  
  c.predictors <- as.matrix(processed$cure_processed$predictors)
  cure_formula <- paste0(" ~ ", paste(colnames(c.predictors), collapse=" + ") )
  cure_formula <- as.formula(cure_formula)
  
  comb.data <- cbind(s.outcomes,s.predictors,c.predictors)
  comb.data <- comb.data[,!duplicated(t(comb.data))]
  newdata <- data.frame(comb.data)
  
  fit <- cureit_impl(surv_formula, cure_formula, newdata, conf.level = conf.level, nboot = nboot, eps = eps)
  
  output <-
    new_cureit(
      surv_coefs = fit$surv_coefs,
      surv_coef_names = fit$surv_coef_names,
      cure_coefs = fit$cure_coefs,
      cure_coef_names = fit$cure_coef_names,
      surv_formula_input = surv_formula_old,
      cure_formula_input = cure_formula_old,
      surv_formula_smcure = surv_formula,
      cure_formula_smcure = cure_formula,
      tidy = fit$tidy,
      smcure = fit$smcure,
      data = data,
      surv_blueprint = processed$surv_processed$blueprint,
      cure_blueprint = processed$cure_processed$blueprint,
      conf.level = conf.level,
      nboot=nboot,
      eps=eps
    )
  
  output
}


# Generic
#' @rdname cureit
#' @export
cureit <- function(object, ...) {
  UseMethod("cureit")
}

# Default
#' @rdname cureit
#' @export
cureit.default <- function(object, ...) {
  stop("`cureit()` is not defined for a '", class(object)[1], "'.", call. = FALSE)
}





