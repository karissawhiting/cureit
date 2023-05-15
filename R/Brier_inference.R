#' Inference for Brier score for cureit objects
#' @rdname Brier_inference_cureit
#' @param object A cureit object.
#' @param times Numeric vector of times to obtain survival probability estimates at
#' @param ... Additional arguments passed to other methods.
#' @return a tibble
#' @family cureit() functions
#' @export
#' @examples

#' trial <- na.omit(trial)
#' p <- cureit(surv_formula = Surv(ttdeath, death) ~ age, 
#'    cure_formula = ~ age,
#'    data = trial) 
#'    
#' times = seq(5,24,0.5)
#' bootbrier <- Brier_inference_bootstrap(p,times=times, nboot = 10)
#' plot(times,bootbrier$brier,type="s",xlab="Time",ylab="Brier score")
#' lines(times,bootbrier$brier_2.5,type="s",lty=2)
#' lines(times,bootbrier$brier_97.5,type="s",lty=2)
#' lines(times,bootbrier$brier_cox,type="s",col="red")
#' lines(times,bootbrier$brier_cox_2.5,type="s",lty=2,col="red")
#' lines(times,bootbrier$brier_cox_97.5,type="s",lty=2,col="red")
#' legend("topleft",c("Cure model","Cox model"),col=c("black","red"),lty=1)

Brier_inference_bootstrap <- function(object, times,...) {
  

  # Data checks -------
  if (object$nboot <= 0) {
    stop("`nboot=` cannot be less than 1. Please rerun `cureit` with `nboot` > 0", call. = FALSE)
  }
  if (!is.null(times) && any(times < 0)) {
    stop("`times=` must be non-negative.", call. = FALSE)
  }

  Brier_original <- predict(
    object = object,
    newdata = object$data,
    method = "prob",
    times = times, 
    brier = TRUE, cox = TRUE
  )


  nboot <- object$nboot
  bootfit <- object$smcure$bootstrap_fit
  nboot_success <- length(bootfit)
  
  # only keep non NULL
  null_index <- purrr::map_lgl(bootfit, is.null)
  
  i <- 1
  
  boot_brier <- boot_brier_cox <- matrix(NA, nrow = nboot_success, ncol = length(times))
  
  # i cannot be larger than min time of any run. maybe make this a try catch?
  while (i <= nboot_success) {
    if(!is.null(bootfit[[i]])) {
      pred_bootfit <- predict(bootfit[[i]],
                              method = "prob",
                              times = times,
                              brier = TRUE,
                              cox = TRUE)

      boot_brier[i, ] <- pred_bootfit$brier
      boot_brier_cox[i, ] <- pred_bootfit$brier_cox
    } else {
      boot_brier[i, ] <- NA
      boot_brier_cox[i, ] <- NA
    }

    i <- i + 1
  }
  
  brier_sd <- apply(boot_brier,2,sd, na.rm = TRUE)
  brier_cox_sd <- apply(boot_brier_cox,2, sd, na.rm = TRUE)
  brier_2.5 <- apply(boot_brier,2,function(x) quantile(x,probs=0.025, na.rm = TRUE))
  brier_97.5 <- apply(boot_brier,2,function(x) quantile(x,probs=0.975, na.rm = TRUE))
  brier_cox_2.5 <- apply(boot_brier_cox,2,function(x) quantile(x,probs=0.025, na.rm = TRUE))
  brier_cox_97.5 <- apply(boot_brier_cox,2,function(x) quantile(x,probs=0.975, na.rm = TRUE))
  
  return(list(brier = Brier_original$brier,
              brier_cox = Brier_original$brier_cox,
              brier_sd = brier_sd,
              brier_cox_sd = brier_cox_sd,
              brier_2.5 = brier_2.5,
              brier_97.5 = brier_97.5,
              brier_cox_2.5 = brier_cox_2.5,
              brier_cox_97.5 = brier_cox_97.5,
              boot_brier = boot_brier,
              boot_brier_cox = boot_brier_cox))
}


#' @rdname Brier_inference_cureit
#' @param object A cureit object.
#' @param ncv Number of cross-validation folds.
#' @param times Numeric vector of times to obtain survival probability estimates at
#' @param ... Additional arguments passed to other methods.
#' @export
#' @examples
#' p <- cureit(surv_formula = Surv(ttdeath, death) ~ age, 
#'    cure_formula = ~ age,
#'    data = trial) 
#'    
#' times = seq(5,24,0.5)
#' cvbrier <- Brier_inference_cv(p,times=times)
#' plot(times,cvbrier$cv_brier_cure_mean,type="s",xlab="Time",ylab="Brier score")
#' lines(times,cvbrier$cv_brier_cox_mean,type="s",lty=2,col="red")

Brier_inference_cv <- function(object,
                               ncv = 10,
                               times,
                               ...) {
  
  # Data checks -------
  if (ncv <= 1){
    stop("`cv=` cannot be less than 2", call. = FALSE)
  }
  if (!is.null(times) && any(times < 0)) {
    stop("`times=` must be non-negative.", call. = FALSE)
  }
  
  processed <- cureit_mold(object$surv_formula, object$cure_formula, object$data,
                           surv_blueprint = object$surv_blueprint, cure_blueprint = object$cure_blueprint)
  s.outcomes <- as.matrix(processed$surv_processed$outcomes[, 1, drop = TRUE])
  status <- s.outcomes[,"status"]
  folds <- createFolds(status, list=FALSE, k = ncv)
  
  brier_cure <- brier_cox <- matrix(NA,nrow=ncv,ncol=length(times))
  
  for (i in 1:ncv){
    
    cv.training <- object$data[folds!=i,]
    cv.testing <- object$data[folds==i,]
    
    cv.fit <- cureit(surv_formula = object$surv_formula,
                     cure_formula = object$cure_formula, 
                     data = cv.training,
                     nboot = 0,
                     eps = 0.05)
    
    cv.predict <- predict(cv.fit,newdata=cv.testing,times=times,brier=TRUE,cox=TRUE)
    brier_cure[i,] <- cv.predict$brier
    brier_cox[i,] <- cv.predict$brier_cox
    
  }
  
  return(list(cv_brier_cure_mean = colMeans(brier_cure),
              cv_brier_cox_mean = colMeans(brier_cox),
              cv_brier_cure_median = apply(brier_cure,2,median),
              cv_brier_cox_median = apply(brier_cox,2,median),
              brier_cure = brier_cure,
              brier_cox = brier_cox))
  
}

