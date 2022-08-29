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
  
  if (survival){
    
    surv.covnames <- all.vars(x$surv_formula[[3]])
    
    
  }
  
  if (cure){
    
    cure.covnames <- all.vars(x$cure_formula[[2]])
    
  }
  
}