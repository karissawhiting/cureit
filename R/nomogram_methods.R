#' Nomogram methods for cureit objects
#'
#' @param survival Logical indicating whether or not to output the nomogram
#' based on the survival submodel. Defaults to `TRUE`.
#' @param cure Logical indicating whether or not to output the nomogram
#' based on the cured probability submodel. Defaults to `TRUE`.
#'
#' @name nomogram_methods_cureit
#' @return a tibble
#' @family cureit() functions
#' @examples
#' cureit <- cureit(Surv(ttdeath, death_cr) ~ age + grade, trial)
#'
#' nomogram(cureit)
NULL

# nomogram
#' @rdname nomogram_methods_cureit
#' @export
#' @family cureit
nomogram.cureit <- function(x,
                            survival = TRUE, 
                            cure = TRUE,
                            ...) {
  
  processed <- cureit_mold(x$surv_formula, x$cure_formula, x$data)
  
  if (survival){
    
    surv.covnames <- all.vars(x$surv_formula[[3]])
    surv.factors <- names(x$surv_xlevels)
    surv.continuous <- setdiff(surv.covnames,surv.factors)
    num.levels <- sapply(x$surv_xlevels,length)
    nmmat <- data.frame(variables=c(unlist(mapply(rep,surv.factors,num.levels)),surv.continuous),
                             levels=c(unlist(x$surv_xlevels),rep(NA,length(surv.continuous))),
                             lp=0,
                             lpscale=0)
    nmmat$levels <- gsub(" ", "_", nmmat$levels)
    nmmat$combined <- ifelse(is.na(nmmat$levels),
                                  nmmat$variables,
                                  paste0(nmmat$variables,nmmat$levels))
    
    covmin <- apply(processed$surv_processed$predictors,2,min)
    covmax <- apply(processed$surv_processed$predictors,2,max)
    lpdiff <- (covmax-covmin)*x$surv_coefs
    nm_scale <- 100/max(abs(lpdiff))
    lpscale <- lpdiff * nm_scale
    
    for (i in 1:length(lpdiff)){
      nmmat$lp[nmmat$combined == names(lpdiff)[i]] = lpdiff[i]
      nmmat$lpscale[nmmat$combined == names(lpdiff)[i]] = lpscale[i]
    }
    
    if (min(nmmat$lp) < 0){
      paste(
        "The effect size must be greater than 0 for categorical variables.",
        "Please adjust the reference group."
      ) %>%
        stop(call. = FALSE)
    }
    
    nmmat$levels <- gsub("_", " ", nmmat$levels)
    
    p1 <- nmmat %>%
      ggplot(aes(x = lpscale, y = variables)) + geom_line() +
      geom_point() + 
      geom_text(aes(label = levels), vjust = 1.5, angle=0)  + ylab(" ") + xlab(" ") +
      ggtitle("Estimated Survival for Uncured") +
      theme_minimal() +
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_blank())

    
  }
  
  if (cure){
    
    cure.covnames <- all.vars(x$cure_formula[[2]])
    
  }
  
}