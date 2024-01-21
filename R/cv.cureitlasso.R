cv.cureitlasso <- function(t,
                           d,
                           x,
                           minmu.ratio=0.05, # minimum penalty value in the logistic model
                           minlambda.ratio=0.05, # minimum penalty value in the cox model
                           length.grid=10,
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
  require(survcomp)
  
  if (progress) print("Fitting by EM algorithm ...")
  
  fit <- cureitlasso(t,d,x,
                     minmu.ratio,
                     minlambda.ratio,
                     length.grid,
                     mus=NULL,
                     lambdas=NULL,
                     tol,
                     maxiter,
                     progress)
  
  if (progress) print("Running cross validations ...")
  
  if (is.null(seed)) seed <- as.numeric(Sys.Date())
  
  # Fold split
  foldid <- coxsplit(as.matrix(Surv(t,d)), nfolds)
  cv_cidx <- array(NA,dim=c(length.grid,length.grid,nfolds))
  
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
          predsurv <- rep(1,length(ti))
          for (l in 1:length(ti)){
            if (t[l] >= min(tj)){
              ids <- which(tj == max(tj[tj <= t[l]]))
              predsurv[l] <- exp(-cumhaz[ids]*predsurvexp[l])
            }
          }
          
          preds <- 1 - predcure + predcure*predsurv
          cv_cidx[j,k,i] <- concordance.index(x=preds,surv.time=ti,surv.event=di)$c.index
          
          
        }
        
      }
      
    }
    
  }else if (ncore > 1){
    
    cl <- makeCluster(ncore)
    registerDoParallel(cl)
    
    cv.fit <- foreach(i = 1:nfolds) %dopar% {
      
      cureitlasso(t[foldid != i],
                  d[foldid != i],
                  x[foldid != i,],
                  minmu.ratio,
                  minlambda.ratio,
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
          predsurv <- rep(1,length(ti))
          for (l in 1:length(ti)){
            if (t[l] >= min(tj)){
              ids <- which(tj == max(tj[tj <= t[l]]))
              predsurv[l] <- exp(-cumhaz[ids]*predsurvexp[l])
            }
          }
          
          preds <- 1 - (1 - predcure + predcure*predsurv)
          cv_cidx[j,k,i] <- concordance.index(x=preds,surv.time=ti,surv.event=di)$c.index
          
        }
        
      }
      
    }
    
  }
  
  cv_cidx_mean <- apply(cv_cidx,c(1,2),function(x) mean(x))
  cv_cidx_se <- apply(cv_cidx,c(1,2),function(x) sd(x)/sqrt(nfolds))
  # pheatmap::pheatmap(cv_cidx_mean,cluster_cols = F,cluster_rows = F)
  
  return(list(fit=fit,
              num_alpha=num_alpha,
              num_beta=num_beta))
  
}