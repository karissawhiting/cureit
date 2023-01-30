

# trial <- na.omit(trial)
# p <- boot_cure(surv_formula = Surv(ttdeath, death) ~ age,
#    cure_formula = ~ age,
#    data = trial, 
#    nboot = 10)
# 
# mod_1 <- p[[1]]$smcure
# predict(p[[1]])







#' Bootstrap smcure fits
#'
#' @param surv_formula 
#' @param cure_formula 
#' @param data 
#' @param conf.level 
#' @param nboot 
#' @param eps 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' trial <- na.omit(trial)
#' p <- boot_cure(surv_formula = Surv(ttdeath, death) ~ age,
#'                cure_formula = ~ age,
#'                data = trial,
#'                nboot = 10)

#' mod_1 <- p[[1]]$smcure
#' 
boot_cure <- function(surv_formula,
                      cure_formula,
                      data,
                      conf.level = 0.95, 
                      nboot = 100,
                      eps = 1e-7,
                      brier,
                      
                    #  object,
                   #   newdata = NULL,
                     # nboot = 100,
                    #  times,
                      ...) {
  
  # Argument checks -----------
  if (nboot <= 0) {
    stop("`nboot=` cannot be less than 1", call. = FALSE)
  }
  
  # Fit cureit object with 0 bootstraps -----
  object <- cureit(surv_formula,
                   cure_formula,
                   data,
                   conf.level = conf.level,
                   nboot = 0,
                   eps = eps, 
                   brier)
  
  # Prep data for bootstrap ----------
  newdata <- object$data
  processed <- cureit_mold(object$surv_formula,
                           object$cure_formula, 
                           object$data,
                           surv_blueprint = object$surv_blueprint, 
                           cure_blueprint = object$cure_blueprint)
  
  s.outcomes <- as.matrix(processed$surv_processed$outcomes[, 1, drop = TRUE])
  status <- s.outcomes[, "status"]
  
  # subset data by endpont status
  data1 <- subset(newdata, status == 1)
  data0 <- subset(newdata, status == 0)
  n1 <- nrow(data1)
  n0 <- nrow(data0)
  
  # Start bootstrapping ----------
  boot_fit_results <- x <- vector("list", nboot)
  
  i <- 1
  
  while (i <= nboot) {
    
    # print(i)
    
    id1 <- sample(1:n1, n1, replace = TRUE)
    id0 <- sample(1:n0, n0, replace = TRUE)
    bootdata <- rbind(data1[id1, ], data0[id0, ])
    
    bootfit <- cureit(
      surv_formula = object$surv_formula,
      cure_formula = object$cure_formula,
      data = bootdata,
      nboot = 0,
      eps = object$eps, 
      brier = brier
    )
    
    boot_fit_results[[i]]  <- bootfit
    print(i)
    i <- i + 1
  }
  
  # returns a list of cureit objects
  boot_fit_results
  
}
  


  