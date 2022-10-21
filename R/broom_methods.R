#' Broom methods for smcure objects
#'
#' @param exponentiate Logical indicating whether or not to exponentiate the
#' coefficient estimates. Defaults to `FALSE`.
#' @param conf.level Level of the confidence interval. Default matches that in
#' `smcure(conf.level=)` (typically, 0.95)
#' @param ... Additional arguments passed to other methods.
#' 
#' @name broom_methods_smcure
#' @return a tibble
#' @family smcure() functions
#' @examples
#' smcure <- smcure::smcure(Surv(ttdeath, death) ~ marker, cureform = ~marker,
#'  data = trial, model = "ph")
#'
#' tidy(smcure)
NULL


# tidy
#' @param x An smcure object created by smcure::smcure()
#' @param conf.int A Logical indicating whether or not to include a confidence interval in the tidied output. Defaults to FALSE.
#' @rdname broom_methods_smcure
#' @export
#' @family smcure tidiers
tidy.smcure <- function(x,
                        exponentiate = FALSE,
                        conf.int = FALSE,
                        conf.level = 0.95, ...) {
  
  # df_logit <- broom::tidy(
  #   x$logistfit,
  #   exponentiate = exponentiate,
  #   conf.int = conf.int,
  #   conf.level = conf.level, ...
  # )
  # 
  # if (isTRUE(conf.int)) {
  #   df_logit <-
  #     df_logit %>%
  #     dplyr::relocate(.data$conf.low, .data$conf.high, .before = .data$p.value)
  # }
  # 
  # df_logit$term <- paste0(x$bnm,", Cure model")
  
  if (conf.int){
    df_logit <- data.frame(term = paste0(x$bnm,", Cure model"),
                           estimate=x$b,
                           std.error=x$b_sd,
                           statistic=x$b_zvalue,
                           conf.low=x$b-1.96*x$b_sd,
                           conf.high=x$b+1.96*x$b_sd,
                           p.value=x$b_pvalue) %>% as_tibble()
    
    df_surv <- data.frame(term=paste0(x$betanm,", Survival model"),
                          estimate=x$beta,
                          std.error=x$beta_sd,
                          statistic=x$beta_zvalue,
                          conf.low=x$beta-1.96*x$beta_sd,
                          conf.high=x$beta+1.96*x$beta_sd,
                          p.value=x$beta_pvalue
    ) %>% as_tibble() 
    if (exponentiate){
      df_logit[,c(2,5,6)] <- exp(df_logit[,c(2,5,6)])
      df_surv[,c(2,5,6)] <- exp(df_surv[,c(2,5,6)])
    }
  }else{
    df_logit <- data.frame(term = paste0(x$bnm,", Cure model"),
                           estimate=x$b,
                           std.error=x$b_sd,
                           statistic=x$b_zvalue,
                           p.value=x$b_pvalue) %>% as_tibble()
    
    df_surv <- data.frame(term=paste0(x$betanm,", Survival model"),
                          estimate=x$beta,
                          std.error=x$beta_sd,
                          statistic=x$beta_zvalue,
                          p.value=x$beta_pvalue) %>% as_tibble() 
    if (exponentiate){
      df_logit[,2] <- exp(df_logit[,2])
      df_surv[,2] <- exp(df_surv[,2])
    }
  }
  
  list(df_cure=df_logit,
       df_surv=df_surv)
  
}

#' Broom methods for cureit objects
#'
#' @param exponentiate Logical indicating whether or not to exponentiate the
#' coefficient estimates. Defaults to `FALSE`.
#' @param conf.level Level of the confidence interval. Default matches that in
#' `smcure(conf.level=)` (typically, 0.95)
#' @inheritParams broom_methods_smcure
#'
#' @name broom_methods_cureit
#' @return a tibble
#' @family cureit() functions
#' @examples
#' cureit <- cureit(surv_formula = Surv(ttdeath, death) ~ marker, cure_formula = ~marker,
#'  data = trial)
#'
#' tidy(cureit)
NULL

# tidy
#' @rdname broom_methods_cureit
#' @export
#' @family cureit tidiers
#' @inheritParams broom_methods_smcure
tidy.cureit <- function(x,
                        exponentiate = FALSE,
                        conf.int = FALSE,
                        conf.level = 0.95, ...) {
  
  # df_logit <- broom::tidy(
  #   x$logistfit,
  #   exponentiate = exponentiate,
  #   conf.int = conf.int,
  #   conf.level = conf.level, ...
  # )
  # 
  # if (isTRUE(conf.int)) {
  #   df_logit <-
  #     df_logit %>%
  #     dplyr::relocate(.data$conf.low, .data$conf.high, .before = .data$p.value)
  # }
  # 
  # df_logit$term <- paste0(x$bnm,", Cure model")
  
  if (conf.int){
    df_logit <- data.frame(term = paste0(x$smcure$bnm,", Cure model"),
                           estimate=x$smcure$b,
                           std.error=x$smcure$b_sd,
                           statistic=x$smcure$b_zvalue,
                           conf.low=x$smcure$b-1.96*x$smcure$b_sd,
                           conf.high=x$smcure$b+1.96*x$smcure$b_sd,
                           p.value=x$smcure$b_pvalue) %>% as_tibble()
    
    df_surv <- data.frame(term=paste0(x$smcure$betanm,", Survival model"),
                          estimate=x$smcure$beta,
                          std.error=x$smcure$beta_sd,
                          statistic=x$smcure$beta_zvalue,
                          conf.low=x$smcure$beta-1.96*x$smcure$beta_sd,
                          conf.high=x$smcure$beta+1.96*x$smcure$beta_sd,
                          p.value=x$smcure$beta_pvalue
    ) %>% as_tibble() 
    if (exponentiate){
      df_logit[,c(2,5,6)] <- exp(df_logit[,c(2,5,6)])
      df_surv[,c(2,5,6)] <- exp(df_surv[,c(2,5,6)])
    }
  }else{
    df_logit <- data.frame(term = paste0(x$smcure$bnm,", Cure model"),
                           estimate=x$smcure$b,
                           std.error=x$smcure$b_sd,
                           statistic=x$smcure$b_zvalue,
                           p.value=x$smcure$b_pvalue) %>% as_tibble()
    
    df_surv <- data.frame(term=paste0(x$smcure$betanm,", Survival model"),
                          estimate=x$smcure$beta,
                          std.error=x$smcure$beta_sd,
                          statistic=x$smcure$beta_zvalue,
                          p.value=x$smcure$beta_pvalue) %>% as_tibble() 
    if (exponentiate){
      df_logit[,2] <- exp(df_logit[,2])
      df_surv[,2] <- exp(df_surv[,2])
    }
  }
  
  # list(df_cure=df_logit,
  #      df_surv=df_surv)
  
  bind_rows(df_logit, df_surv)
  
}


