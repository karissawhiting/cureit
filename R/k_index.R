
#' K-index estimation
#'
#' @export
#' @param object a cureit model object or a coxph model object
.calc_k_index <- function(object){
  
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
                             object$data,
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
          
          N <- sum(logit_inv_cure)
          comp_risk_cox <- outer(lp_surv_model, lp_surv_model, .kernel_fun)
          diag(comp_risk_cox) <- 0
          
          # Based on definition
          sum(comp_risk_cox * outer(logit_inv_cure, logit_inv_cure, "*")  ) * 2  / ( N^2 - sum(logit_inv_cure^2)) 

         },
         "coxph" = {
           
           lp_surv_model = predict(object, type = "risk") 
           # set cure probablity always be 1
           logit_inv_cure  <- rep(1, object$n)        
           
           lp_surv_model <- as.numeric(lp_surv_model)
           logit_inv_cure <- as.numeric(logit_inv_cure)
           
           N <- sum(logit_inv_cure)
           comp_risk_cox <- outer(lp_surv_model, lp_surv_model, .kernel_fun)
           diag(comp_risk_cox) <- 0
           
           # Based on definition
           sum(comp_risk_cox * outer(logit_inv_cure, logit_inv_cure, "*")  ) * 2  / ( N^2 - sum(logit_inv_cure^2)) 
         })
  k
}


#' K-index estimation
#'
#' @export
#' @param object a cureit model object
k_index <- function(object){
  
# Bootstrap CIs -----------------------------------------------------------

  nboot <- object$nboot
  
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
  boot_kindex_2.5 <- apply(boot_kindex,2,function(x) quantile(x,probs=0.025))
  boot_kindex_97.5 <- apply(boot_kindex,2,function(x) quantile(x,probs=0.975))
  boot_kindex_mean <- apply(boot_kindex,2,function(x) mean(x, na.rm = TRUE))

  # compare to cox model ------------------
  if (compare_cox){
    cox_marginal = probs_at_times(scox_prd,times)
    
    brier_cox <- colSums(mapply(function(count,ipw,surv_marginal,obs) count*(obs-surv_marginal)^2/(ipw+1e-4),
                                count,ipw,cox_marginal,obs))/nrow(s.outcomes)
    
    names(brier_cox) <- paste0("T = ",times)
    
    return(list(cured = cure_prd,
                surv_uncured = probs_at_times(scure_prd,times),
                surv_marginal = probs_at_times(spop_prd,times),
                brier = brier,
                brier_cox = brier_cox
    ))
    
  return(list(k_index = k,
              boot_sd = boot_kindex_sd,
              boot_kindex_2.5 = boot_kindex_2.5,
              boot_kindex_97.5 = boot_kindex_97.5,
              boot_kindex_mean = boot_kindex_mean,
              boot_kindex = boot_kindex))  
  } else {
    .calc_k_index(object)
  }

}



