times = seq(5,24,0.5)

#' object <- cureit(surv_formula = Surv(ttdeath, death) ~ age, 
#'    cure_formula = ~ age,
#'    data = trial) 
#'    

brier <- function(object, times = NULL, probs = NULL,
                           newdata = NULL, method="prob",
                           brier = FALSE, cox = FALSE, ...) {
  

  # Get Necessary Values ----------------------------------------------------

  
  # getting predictions on the original model fit ------------------------------
  
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
  
  if (is.vector(newX)) {
    newX = as.matrix(newX)
  }
  
  s0 = as.matrix(object$smcure$s[order(object$smcure$Time)], ncol = 1)
  n = nrow(s0)
  uncureprob = exp(object$smcure$b %*% t(newZ))/(1 + exp(object$smcure$b %*% t(newZ)))
  
  scure = array(0, dim = c(n, nrow(newX)))
  t = array(0, dim = c(n, nrow(newX)))
  spop = array(0, dim = c(n, nrow(newX)))
  
  #### Only support PH model for now
  ebetaX = exp(object$smcure$beta %*% t(newX))
  for (i in 1:nrow(newZ)) {
    scure[, i] = s0^ebetaX[i]
  }
  
  for (i in 1:n) {
    for (j in 1:nrow(newX)) {
      spop[i, j] = uncureprob[j] * scure[i, j] + (1 - 
                                                    uncureprob[j])
    }
  }
  
  spop_prd = cbind(object$smcure$Time, spop)
  scure_prd = cbind(object$smcure$Time, scure)
  cure_prd = 1-uncureprob
  
  # Brier -------------------------------------------------------------------

  
    s.outcomes <- as.matrix(processed$surv_processed$outcomes[, 1, drop = TRUE])
    s.outcomes[,"status"] <- 1 - s.outcomes[,"status"] # censoring indicator
    cens_formula <- paste0("Surv(", paste(colnames(s.outcomes),collapse=","), ")"," ~ 1")
    cens_formula <- as.formula(cens_formula)
    cens_fit <- survfit(cens_formula,data.frame(s.outcomes))
    cens_prd <- cbind(cens_fit$time,cens_fit$surv)
    s.outcomes[,"status"] <- 1 - s.outcomes[,"status"] # transform it back to event indicator
    
    if (cox){
      newdat <- newdata %||% object$data
      coxdat <<- object$data
      coxformula <- object$surv_formula
      coxfit <- coxph(coxformula,coxdat)
      s0 <- survfit(coxfit)$surv
      scox = array(0, dim = c(length(s0), nrow(newdat)))
      coxlp <- predict(coxfit,newdata=newdat,type="lp")
      for (i in 1:nrow(newdat)){
        scox[,i] <- s0^exp(coxlp[i])
      }
      scox_prd <- cbind(survfit(coxfit)$time,scox)
    }
    


    # Must be Probabilities ----
    surv_marginal = probs_at_times(spop_prd, times)
    na <- purrr::map(surv_marginal, ~which(is.na(.x)))
    na_index <- na[[1]]
    
    surv_marginal <- lapply(surv_marginal, function(x) {
      
      if(length(na_index) > 0 ) {
        x <- x[-na_index]
      }
      
      x
    })
    
    count <- lapply(times, function(x) { 
      r <- ifelse(s.outcomes[,"time"] > x, 1, s.outcomes[,"status"])
      if(length(na_index) > 0 ) {
        r <- r[-na_index]
      }
      r
    })
    
    ipw <- lapply(times,function(x) {
      r <- ifelse(s.outcomes[,"time"] > x, unlist(probs_at_times(cens_prd,x)), cens_prd[,2])
      if(length(na_index) > 0 ) {
        r <- r[-na_index]
      }
      r
    })
    
    obs <- lapply(times, function(x) {
      r <- ifelse(s.outcomes[,"time"] > x, 1, 0)
      if(length(na_index) > 0 ) {
        r <- r[-na_index]
      }
      r
    })
    
    brier <- colSums(mapply(function(count,ipw,surv_marginal,obs) count*(obs-surv_marginal)^2/(ipw+1e-4),
                            count,ipw,surv_marginal,obs), na.rm = TRUE)/nrow(s.outcomes)
    names(brier) <- paste0("T = ",times)
    
    if (cox){
      cox_marginal = probs_at_times(scox_prd,times)
      
      brier_cox <- colSums(mapply(function(count,ipw,surv_marginal,obs) count*(obs-surv_marginal)^2/(ipw+1e-4),
                                  count,ipw,cox_marginal,obs), na.rm = TRUE)/nrow(s.outcomes)
      names(brier_cox) <- paste0("T = ",times)
      
      ret <- list(cured = cure_prd,
                  surv_uncured = probs_at_times(scure_prd,times),
                  surv_marginal = probs_at_times(spop_prd,times),
                  brier = brier,
                  brier_cox = brier_cox
      )
      
    }else{
      
      ret <- list(cured = cure_prd,
                  surv_uncured = probs_at_times(scure_prd,times),
                  surv_marginal = probs_at_times(spop_prd,times),
                  brier = brier
      )
      
    }
    
}




brier_survival_impl <- function(truth,
                                estimate,
                                censoring_weights,
                                case_weights,
                                eval_time) {
  
  surv_time <- .extract_surv_time(truth)
  surv_status <- .extract_surv_status(truth)
  
  if (!is.null(case_weights)) {
    norm_const <- sum(case_weights)
  } else {
    case_weights <- rep(1, length(estimate))
    norm_const <- sum(!survival::is.na.Surv(truth))
  }
  
  category_1 <- surv_time <= eval_time & surv_status == 1
  category_2 <- surv_time > eval_time
  
  # (0 - estimate) ^ 2 == estimate ^ 2
  res <- (category_1 * estimate^2 * censoring_weights) +
    (category_2 * (1 - estimate)^2 * censoring_weights)
  
  res <- res * case_weights
  res <- sum(res, na.rm = TRUE)
  res / norm_const
}
  


