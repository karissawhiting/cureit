
#' K-index estimation
#'
#' @export
#' @param object a cureit model object or a coxph model object
#' @param newdata A `base::data.frame()` or `tibble::tibble()` containing all
#' the original predictors used to create object. Defaults to `NULL`.
#' @examples 
#' 
#' p <- cureit(surv_formula = Surv(ttdeath, death) ~ age, 
#'    cure_formula = ~ age,
#'    data = trial) 
#' .calc_k_index(p)
#' 
#' cox_p <- survival::coxph(Surv(ttdeath, death) ~ age,
#'    data = trial)
#'
#' # Need to debug this
#' # .calc_k_index(cox_p)
#' 
.calc_k_index <- function(object, newdata = NULL){
  
  model_class <- if(inherits(object, "cureit")) {
    "cureit"
  } else if(inherits(object, "coxph")) {
    "coxph"
  }
  
  model_class %||% 
    cli::cli_abort("Model must be an object of type `cureit` or `coxph`")
  
  k = switch(model_class, 
             "cureit" = {
               processed <- cureit_mold(object$surv_formula,
                                        object$cure_formula,
                                        newdata %||% object$data,
                                        surv_blueprint = object$surv_blueprint,
                                        cure_blueprint = object$cure_blueprint)
               
               newX=processed$surv_processed$predictors
               newZ=processed$cure_processed$predictors
               
               if (is.vector(newZ)) {
                 newZ = as.matrix(newZ)
               }
               newZ = cbind(1, newZ)
               
               if (is.vector(newX))  {
                 newX = as.matrix(newX)
               }
               
               lp_surv_model = t(object$smcure$beta %*% t(newX))
               lp_cure_model = t(object$smcure$b %*% t(newZ))
               logit_inv_cure  <- .logit_inv(lp_cure_model)
               
               lp_surv_model <- as.numeric(lp_surv_model)
               logit_inv_cure <- as.numeric(logit_inv_cure)
               
               N <- sum(logit_inv_cure, na.rm = TRUE)
               comp_risk_cox <- outer(lp_surv_model, lp_surv_model, .kernel_fun)
               diag(comp_risk_cox) <- 0
               
               # Based on definition
               sum(comp_risk_cox * outer(logit_inv_cure, logit_inv_cure, "*"), na.rm = TRUE  ) * 2  / ( N^2 - sum(logit_inv_cure^2, na.rm = TRUE)) 
               
             },
             "coxph" = {
               
               lp_surv_model = predict(object, type = "risk", newdata=newdata %||% object$data) 
               # set cure probablity to always be 1 in the case of a coxmodel (source: evacure package)
               logit_inv_cure  <- rep(1, length(lp_surv_model)) # this would address the sample size discrepancy when newdata is supplied        
               
               lp_surv_model <- as.numeric(lp_surv_model)
               logit_inv_cure <- as.numeric(logit_inv_cure)
               
               N <- sum(logit_inv_cure)
               comp_risk_cox <- outer(lp_surv_model, lp_surv_model, .kernel_fun)
               diag(comp_risk_cox) <- 0
               
               # Based on definition
               sum(comp_risk_cox * outer(logit_inv_cure, logit_inv_cure, "*")  ) * 2  / ( N^2 - sum(logit_inv_cure^2)) 
             }
             )

  k
}


#' K-index estimation
#'
#' @export
#' @param object a cureit model object
#' @param newdata A `base::data.frame()` or `tibble::tibble()` containing all
#' the original predictors used to create object. Defaults to `NULL`.
#' @param ncv number of folds for cross validation. Defaults to `NULL`. 
#' If specified, `newdata` will be ignored and cross-validated k-index will be reported.
#' @examples 
#' p <- cureit(surv_formula = Surv(ttdeath, death) ~ age, 
#'    cure_formula = ~ age,
#'    data = trial) 
#' k_index(p)
k_index <- function(object, newdata=NULL, ncv=NULL){
  
  # Bootstrap CIs -----------------------------------------------------------
  
  k = .calc_k_index(object, newdata)
  nboot <- object$nboot
  
  if(!is.null(ncv)){
    
    # Data checks -------
    if (ncv <= 1){
      stop("`cv=` cannot be less than 2", call. = FALSE)
    }
    
    if(nboot > 1) {
      bootfit <- object$smcure$bootstrap_fit
      i <- 1
      boot_kindex <- matrix(NA, nrow = nboot, ncol = 1)
      
      while (i <= nboot) {
        pred_bootfit <- .calc_k_index(bootfit[[i]])
        boot_kindex[i, ] <- pred_bootfit
        
        i <- i + 1
      }
      
      boot_kindex_sd <- apply(boot_kindex, 2, sd)
      boot_kindex_2.5 <- apply(boot_kindex,2,function(x) quantile(x,probs=0.025, na.rm = TRUE))
      boot_kindex_97.5 <- apply(boot_kindex,2,function(x) quantile(x,probs=0.975, na.rm = TRUE))
      boot_kindex_mean <- apply(boot_kindex,2,function(x) mean(x, na.rm = TRUE))
    }
    
    ### Cross validation
    
    processed <- cureit_mold(object$surv_formula, object$cure_formula, object$data,
                             surv_blueprint = object$surv_blueprint, cure_blueprint = object$cure_blueprint)
    s.outcomes <- as.matrix(processed$surv_processed$outcomes[, 1, drop = TRUE])
    status <- s.outcomes[,"status"]
    folds <- createFolds(status, list=FALSE, k = ncv)
    
    cv_kindex <- matrix(NA, nrow = ncv, ncol = 1)
    cv_cox_kindex <- matrix(NA, nrow = ncv, ncol = 1)
    
    for (i in 1:ncv){
      
      cv.training <- object$data[folds!=i,]
      cv.testing <- object$data[folds==i,]
      
      cv.fit <- cureit(surv_formula = object$surv_formula,
                       cure_formula = object$cure_formula, 
                       data = cv.training,
                       nboot = 0,
                       eps = object$eps)
      cv_kindex[i,] <- .calc_k_index(cv.fit, newdata=cv.testing)
      
      cv.coxfit <- coxph(object$surv_formula,data=cv.training)
      cv_cox_kindex[i,] <- .calc_k_index(cv.coxfit, newdata=cv.testing)
      
    }
    
    return(list(k_index = k,
                boot_sd = boot_kindex_sd,
                boot_kindex_2.5 = boot_kindex_2.5,
                boot_kindex_97.5 = boot_kindex_97.5,
                boot_kindex_mean = boot_kindex_mean,
                boot_kindex = boot_kindex,
                cv_kindex = cv_kindex,
                cv_cox_kindex = cv_cox_kindex))  
  }else{
    
    if(is.null(newdata)){
      if(nboot > 1) {
        bootfit <- object$smcure$bootstrap_fit
        i <- 1
        boot_kindex <- matrix(NA, nrow = nboot, ncol = 1)
        
        while (i <= nboot) {
          pred_bootfit <- .calc_k_index(bootfit[[i]])
          boot_kindex[i, ] <- pred_bootfit
          
          i <- i + 1
        }
        
        boot_kindex_sd <- apply(boot_kindex, 2, sd)
        boot_kindex_2.5 <- apply(boot_kindex,2,function(x) quantile(x,probs=0.025, na.rm = TRUE))
        boot_kindex_97.5 <- apply(boot_kindex,2,function(x) quantile(x,probs=0.975, na.rm = TRUE))
        boot_kindex_mean <- apply(boot_kindex,2,function(x) mean(x, na.rm = TRUE))
        
        return(list(k_index = k,
                    boot_sd = boot_kindex_sd,
                    boot_kindex_2.5 = boot_kindex_2.5,
                    boot_kindex_97.5 = boot_kindex_97.5,
                    boot_kindex_mean = boot_kindex_mean,
                    boot_kindex = boot_kindex))  
      }
      
    }else {
      return(k)
    }
  }
  
  
}




