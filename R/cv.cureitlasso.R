# Returns:
# fit - 
# arg: length.grid - 

cv.cureitlasso <- function(t,
                           d,
                           x,
                           minmu.ratio=0.05, # minimum penalty value in the logistic model
                           minlambda.ratio=0.05, # minimum penalty value in the cox model
                           adaptive=FALSE,
                           length.grid=10, # tells how many hyperparameters to test. default is 10 for each 
                           nfolds=5,
                           tol=1e-2,
                           maxiter=100,
                           ncore=1,
                           seed=NULL,
                           progress=FALSE
){
  
  require(glmnet)
  require(survival)
  require(foreach)
  require(parallel)
  require(doParallel)
  require(doSNOW)
  require(survcomp)
  require(pracma)
  
  # order data by times 
  order_idx <- order(t)
  t <- t[order_idx]
  d <- d[order_idx]
  x <- x[order_idx,]
  
  if (progress) print("Fitting by EM algorithm ...")
  
  # main model with no cross validation
  fit <- cureitlasso(t,d,x,
                     minmu.ratio,
                     minlambda.ratio,
                     adaptive,
                     length.grid,
                     mus=NULL,
                     lambdas=NULL,
                     tol,
                     maxiter,
                     progress = TRUE)
  
  if (progress) print("Running cross validations ...")
  
  if (is.null(seed)) seed <- as.numeric(Sys.Date())
  
  # Fold split
  foldid <- coxsplit(as.matrix(Surv(t,d)), nfolds)
  
  # grid/array of values 
  cv_brier <- array(NA,dim=c(length.grid,length.grid,nfolds))
  
  # Run CVs
  if (ncore == 1){
    
    cv.fit <- list()
    
    for (i in 1:nfolds){
      
      if (progress) print(i)
      
      
      cv.fit[[i]] <- cureitlasso(t[foldid != i],
                                 d[foldid != i],
                                 x[foldid != i,],
                                 minmu.ratio,
                                 minlambda.ratio,
                                 adaptive,
                                 length.grid,
                                 mus=fit$mus,
                                 lambdas=fit$lambdas,
                                 tol,
                                 maxiter,
                                 progress)
      
      # Model assessment for CVs
      for (j in 1:length.grid){
        
        for (k in 1:length.grid){
          
          ti <- t[foldid==i]
          di <- d[foldid==i]
          tj <- cv.fit[[i]]$fit[[j]][[k]]$tj
          
          fitcure <- cv.fit[[i]]$fit[[j]][[k]]$fitcure
          predcure <- predict(fitcure,newx=x[foldid==i,],s=min(fitcure$lambda),type="response")
          fitcox <- cv.fit[[i]]$fit[[j]][[k]]$fitcox
          predsurvexp <- predict(fitcox,newx=x[foldid==i,],s=min(fitcox$lambda),type="response")
          haz <- cv.fit[[i]]$fit[[j]][[k]]$haz
          cumhaz <- cv.fit[[i]]$fit[[j]][[k]]$cumhaz
          
          predcens <- survfit(Surv(ti,1-di)~1)$surv
          tcens <- survfit(Surv(ti,1-di)~1)$time
          
          tbrier <- sort(ti)
          brier <- rep(NA,length(tbrier))
          
          for (l in 1:length(tbrier)){
            
            predsurv <- rep(NA,length(ti)) # Conditional survival prob for all patients at tbrier[l]
            
            if (tbrier[l] >= min(tj)){
              
              ids <- which(tj == max(tj[tj <= tbrier[l]]))
              predsurv <- exp(-cumhaz[ids]*predsurvexp)
              ipw <- as.numeric(ti > tbrier[l])*predcens[l] + as.numeric(ti <= tbrier[l])*predcens
              
            }else if (tbrier[l] < min(tj)){
              
              predsurv <- rep(1,length(ti))
              ipw <- 1
              
            }
            
            preds <- 1 - predcure + predcure*predsurv
            brier[l] <- mean(I(di==0)*(1 - preds)^2/(ipw+0.001) + I(di==1)*(0 - preds)^2/(ipw+0.001))
            
          }
          
          # calcualte trapazoidal area for brier score 
          cv_brier[j,k,i] <- trapz(tbrier,brier)
          
          
        }
        
      }
      
    }
    
  }else if (ncore > 1){
    
    cl <- makeCluster(ncore)
    registerDoParallel(cl)
    
    cv.fit <- foreach(i = 1:nfolds) %dopar% {
      
      source(here::here("R/cureitlasso.R"))
      
      cureitlasso(t[foldid != i],
                  d[foldid != i],
                  x[foldid != i,],
                  minmu.ratio,
                  minlambda.ratio,
                  adaptive,
                  length.grid,
                  mus=fit$mus,
                  lambdas=fit$lambdas,
                  tol,
                  maxiter)
      
    }
    
    stopCluster(cl)
    
    # Model assessment for CVs
    
    for (i in 1:nfolds){
      
      for (j in 1:length.grid){
        
        for (k in 1:length.grid){
          
          ti <- t[foldid==i]
          di <- d[foldid==i]
          tj <- cv.fit[[i]]$fit[[j]][[k]]$tj
          
          fitcure <- cv.fit[[i]]$fit[[j]][[k]]$fitcure
          predcure <- predict(fitcure,newx=x[foldid==i,],s=min(fitcure$lambda),type="response")
          fitcox <- cv.fit[[i]]$fit[[j]][[k]]$fitcox
          predsurvexp <- predict(fitcox,newx=x[foldid==i,],s=min(fitcox$lambda),type="response")
          haz <- cv.fit[[i]]$fit[[j]][[k]]$haz
          cumhaz <- cv.fit[[i]]$fit[[j]][[k]]$cumhaz
          
          predcens <- survfit(Surv(ti,1-di)~1)$surv
          tcens <- survfit(Surv(ti,1-di)~1)$time
          
          tbrier <- sort(ti)
          brier <- rep(NA,length(tbrier))
          
          for (l in 1:length(tbrier)){
            
            predsurv <- rep(NA,length(ti)) # Conditional survival prob for all patients at tbrier[l]
            
            if (tbrier[l] >= min(tj)){
              
              ids <- which(tj == max(tj[tj <= tbrier[l]]))
              predsurv <- exp(-cumhaz[ids]*predsurvexp)
              ipw <- as.numeric(ti > tbrier[l])*predcens[l] + as.numeric(ti <= tbrier[l])*predcens
              
            }else if (tbrier[l] < min(tj)){
              
              predsurv <- rep(1,length(ti))
              ipw <- 1
              
            }
            
            preds <- 1 - predcure + predcure*predsurv
            brier[l] <- mean(I(ti > tbrier[l])*(1 - preds)^2/(ipw+0.001) + I(di==1 & ti <= tbrier[l])*(0 - preds)^2/(ipw+0.001))
            
          }
          
          cv_brier[j,k,i] <- trapz(tbrier[1:(l-1)],brier[1:(l-1)])
          
          
        }
        
      }
      
    }
    
  }
  
  # summarize - take average across all folds 
  cv_brier_mean <- apply(cv_brier,c(1,2),function(x) mean(x))
  cv_brier_se <- apply(cv_brier,c(1,2),function(x) sd(x)/sqrt(nfolds))
  # pheatmap::pheatmap(cv_brier_mean,cluster_cols = F,cluster_rows = F)
  
  idxmin <- which(cv_brier_mean == min(cv_brier_mean), arr.ind = TRUE)
  brier1se <- suppressWarnings(min(cv_brier_mean[1:idxmin[1],1:idxmin[2]][which(cv_brier_mean[1:idxmin[1],1:idxmin[2]] > cv_brier_mean[idxmin]+cv_brier_se[idxmin])]))
  if (!is.infinite(brier1se)){
    idx1se <- which(cv_brier_mean==brier1se, arr.ind = TRUE)
  }else{
    brier1se <- max(cv_brier_mean[1:idxmin[1],1:idxmin[2]])
    idx1se <- which(cv_brier_mean==brier1se, arr.ind = TRUE)
  }
  
  selectedmin <- list(cure=which(fit$fit[[idxmin[1]]][[idxmin[2]]]$alpha!=0),
                      cox=which(fit$fit[[idxmin[1]]][[idxmin[2]]]$beta!=0))
  
  selected1se <- list(cure=which(fit$fit[[idx1se[1]]][[idx1se[2]]]$alpha!=0),
                      cox=which(fit$fit[[idx1se[1]]][[idx1se[2]]]$beta!=0))
  
  return(list(fit=fit$fit,
              mus=fit$mus,
              lambdas=fit$lambdas,
              cv_brier_mean = cv_brier_mean,
              cv_brier_se = cv_brier_se,
              index = list(min=idxmin,
                           `1se`= idx1se), 
                           `1se`= idx1se),
              foldid = foldid
  )
  
  
}
