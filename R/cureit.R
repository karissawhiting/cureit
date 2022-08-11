#' Cure Model Regression
#'
#' @param x input object
#' @param formula formula with `Surv()` on LHS and covariates on RHS.
#' @param data data frame
#' @param conf.level confidence level. Default is 0.95.
#' @param nboot number of bootstrap samples used for inference.
#' @param ... passed to methods
#'
#' @return cureit object
#' @family cureit() functions
#' @name cureit
#' @examples
#' cureit(Surv(ttdeath, death_cr) ~ age + grade, trial)
NULL

# Formula method
#' @rdname cureit
#' @export
cureit.formula <- function(formula, data, conf.level = 0.95, nboot = 100, ...) {
  
  # checking inputs  -------------------------
  checking(formula = formula, data = data)
  
  # process model variables ----------------------------------------------------
  processed <- cureit_mold(formula, data)
  
  # building model -------------------------------------------------------------
  cureit_bridge(processed, formula, data, conf.level = conf.level, nboot = nboot)
}

cureit_mold <- function(formula, data) {
  processed <-
    hardhat::mold(
      formula, data,
      blueprint = hardhat::default_formula_blueprint(intercept = TRUE)
    )
  # remove intercept
  processed$predictors <- processed$predictors[, -1]
  processed$predictors <- processed$predictors %>% janitor::clean_names()
  processed$outcomes <- processed$outcomes %>% janitor::clean_names()
  
  processed
}

checking <- function(formula, data, keep_all = FALSE) {
  # evaluating LHS of formula --------------------------------------------------
  formula_lhs <-
    tryCatch(
      {
        rlang::f_lhs(formula) %>%
          rlang::eval_tidy(data = data)
      },
      error = function(e) {
        cli::cli_alert_danger("There was an error evaluating the LHS of the formula.")
        stop(e, call. = FALSE)
      }
    )
  
  # checking type of LHS -------------------------------------------------------
  if (!inherits(formula_lhs, "Surv") ||
      !identical(attr(formula_lhs, "type"), "right")) {
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

new_cureit <- function(coefs, coef_names, surv_formula, cure_formula, tidy, smcure, data,
                    blueprint, conf.level, nboot) {
  
  # function to create an object
  
  if (!is.numeric(coefs)) {
    stop("`coefs` should be a numeric vector.", call. = FALSE)
  }
  
  if (!is.character(coef_names)) {
    stop("`coef_names` should be a character vector.", call. = FALSE)
  }
  
  if (length(coefs) != length(coef_names)) {
    stop("`coefs` and `coef_names` must have the same length.")
  }
  
  hardhat::new_model(
    coefs = coefs %>% stats::setNames(coef_names),
    surv_formula = surv_formula,
    data = data,
    conf.level = conf.level,
    nboot = nboot,
    xlevels =
      stats::model.frame(surv_formula, data = data)[, -1, drop = FALSE] %>%
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
    blueprint = blueprint,
    class = "cureit"
  )
}

cureit_impl <- function(surv_formula, cure_formula, newdata, conf.level = conf.level, nboot = nboot) {
  
  # function to run cureit and summarize with tidy (implementation)
  cureit_fit <-
    smcure::smcure(formula=surv_formula,
                   cureform=cure_formula,
                   data=newdata,
                   model="ph",
                   nboot=nboot
    )
  
  # broom method can be constructed later 
  tidy <- broom::tidy(cureit_fit, conf.int = TRUE, conf.level = conf.level)
  
  coefs <- tidy$estimate
  coef_names <- tidy$term
  
  list(
    coefs = coefs,
    coef_names = coef_names,
    tidy = tidy,
    smcure = cureit_fit
  )
}

cureit_bridge <- function(processed, formula, data, conf.level, nboot) {
  
  # function to connect object and implementation
  
  # validate_outcomes_are_univariate(processed$outcomes)
  
  predictors <- as.matrix(processed$predictors)
  outcomes <- as.matrix(processed$outcomes[, 1, drop = TRUE])
  newdata <- data.frame(cbind(outcomes,predictors))
  new_formula <- paste0("Surv(", paste(colnames(outcomes),collapse=","), ")",
                        " ~ ", paste(colnames(predictors), collapse=" + ") )
  surv_formula <- as.formula(new_formula)
  # NOTE: eventually need to allow different cure variables 
  cure_formula <- surv_formula[-2]
  
  fit <- cureit_impl(surv_formula, cure_formula, newdata, conf.level = conf.level, nboot = nboot)
  
  output <-
    new_cureit(
      coefs = fit$coefs,
      coef_names = fit$coef_names,
      surv_formula = surv_formula,
      cure_formula = cure_formula,
      data = newdata,
      tidy = fit$tidy,
      smcure = fit$smcure,
      blueprint = processed$blueprint,
      conf.level = conf.level,
      nboot=nboot
    )
  
  output
}


# Generic
#' @rdname cureit
#' @export
cureit <- function(x, ...) {
  UseMethod("cureit")
}

# Default
#' @rdname cureit
#' @export
cureit.default <- function(x, ...) {
  stop("`cureit()` is not defined for a '", class(x)[1], "'.", call. = FALSE)
}