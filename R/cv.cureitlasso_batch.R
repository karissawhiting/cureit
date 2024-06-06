cv.cureitlasso_batch <- function(t,
                                 d,
                                 x,
                                 b,
                                 minmu.ratio=0.1, # minimum penalty value in the logistic model
                                 minlambda.ratio=0.1, # minimum penalty value in the cox model
                                 adaptive=FALSE,
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
  require(doSNOW)
  require(survcomp)
  require(pracma)
  
  order_idx <- order(t)
  t <- t[order_idx]
  d <- d[order_idx]
  x <- x[order_idx,]
  b <- b[order_idx]
  
  if (progress) print("Fitting by EM algorithm ...")
  
  fit <- cureitlasso_batch(t,d,x,b,
                           minmu.ratio,
                           minlambda.ratio,
                           adaptive,
                           length.grid,
                           mus=NULL,
                           lambdas=NULL,
                           tol,
                           maxiter,
                           progress)
  
  if (progress) print("Running cross validations ...")
  
  if (is.null(seed)) seed <- as.numeric(Sys.Date())
  
  # Fold split
  foldid <- coxsplitb(as.matrix(Surv(t,d)), b, nfolds)
  cv_brier <- array(NA,dim=c(length.grid,length.grid,nfolds))
  
  # Run CVs
  if (ncore == 1){
    
    cv.fit <- list()
    
    for (i in 1:nfolds){
      
      if (progress) print(i)
      
      cv.fit[[i]] <- cureitlasso_batch(t[foldid != i],
                                       d[foldid != i],
                                       x[foldid != i,],
                                       b[foldid != i],
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
          bi <- b[foldid==i]
          
          num.batches <- length(unique(bi))
          
          vec.d <- rep(NA,nrow(x[foldid==i,]))
          for (ii in 1:num.batches){
            vec.d[bi==ii] <- sum(di[bi == ii])
          }
          
          fitcure <- cv.fit[[i]]$fit[[j]][[k]]$fitcure
          #predcure <- predict(fitcure,newx=x[foldid==i,],s=min(fitcure$lambda),type="response")
          
          product <- exp(x[foldid==i,] %*% fitcure$beta[length(fitcure$lambda),]) # rows: samples; columns: penalty parameters
          predcure <- product
          for (kk in 1:num.batches){
            predcure[bi==kk] <- product[bi==kk]*vec.d[bi==kk]/sum(product[bi==kk])
            predcure[bi==kk] <- predcure[bi==kk] * exp(-predcure[bi==kk])
          }
          
          fitcox <- cv.fit[[i]]$fit[[j]][[k]]$fitcox
          predsurvexp <- predict(fitcox,newx=x[foldid==i,],s=min(fitcox$lambda),type="response")
          haz <- cv.fit[[i]]$fit[[j]][[k]]$haz
          cumhaz <- cv.fit[[i]]$fit[[j]][[k]]$cumhaz
          
          predcens <- survfit(Surv(ti,1-di)~1)$surv
          # tcens <- survfit(Surv(ti,1-di)~1)$time
          
          tbrier <- sort(ti)
          brier <- rep(NA,length(tbrier))
          
          for (l in 1:length(tbrier)){
            
            predsurv <- ipw <- rep(NA,length(ti)) # Conditional survival prob for all patients at tbrier[l]
            
            for (ii in 1:num.batches){
              
              bidx <- which(bi == ii)
              
              if (tbrier[l] < min(tj[[ii]])){
                
                predsurv[bidx] <- 1
                ipw[bidx] <- 1
                
              }else if (tbrier[l] >= min(tj[[ii]])){
                
                ids <- which(tj[[ii]] == max(tj[[ii]][tj[[ii]] <= tbrier[l]]))
                predsurv[bidx] <- exp(-cumhaz[[ii]][ids]*predsurvexp[bidx])
                ipw[bidx] <- as.numeric(ti[bidx] > tbrier[l])*predcens[l] + as.numeric(ti[bidx] <= tbrier[l])*predcens[bidx]
                
              }
              
            }
            
            preds <- 1 - predcure + predcure*predsurv
            brier[l] <- mean(I(ti > tbrier[l])*(1 - preds)^2/(ipw+0.001) + I(di==1 & ti <= tbrier[l])*(0 - preds)^2/(ipw+0.001))
            
          }
          
          # plot(tbrier,brier,type="S",ylim=c(0,1))
          
          cv_brier[j,k,i] <- trapz(tbrier,brier)
          
          
        }
        
      }
      
    }
    
  }else if (ncore > 1){
    
    cl <- makeCluster(ncore)
    registerDoParallel(cl)
    
    cv.fit <- foreach(i = 1:nfolds) %dopar% {
      
      source("~/Projects/Whiting-Qin-cureit/cureit/R/cureitlasso_batch.R")
      
      cureitlasso_batch(t[foldid != i],
                        d[foldid != i],
                        x[foldid != i,],
                        b[foldid != i],
                        minmu.ratio,
                        minlambda.ratio,
                        adaptive,
                        length.grid,
                        mus=fit$mus,
                        lambdas=fit$lambdas,
                        tol,
                        maxiter,
                        progress)
      
    }
    
    stopCluster(cl)
    
    # Model assessment for CVs
    
    for (i in 1:nfolds){
      
      for (j in 1:length.grid){
        
        for (k in 1:length.grid){
          
          ti <- t[foldid==i]
          di <- d[foldid==i]
          tj <- cv.fit[[i]]$fit[[j]][[k]]$tj
          bi <- b[foldid==i]
          
          num.batches <- length(unique(bi))
          
          vec.d <- rep(NA,nrow(x[foldid==i,]))
          for (ii in 1:num.batches){
            vec.d[bi==ii] <- sum(di[bi == ii])
          }
          
          fitcure <- cv.fit[[i]]$fit[[j]][[k]]$fitcure
          #predcure <- predict(fitcure,newx=x[foldid==i,],s=min(fitcure$lambda),type="response")
          
          product <- exp(x[foldid==i,] %*% fitcure$beta[length(fitcure$lambda),]) # rows: samples; columns: penalty parameters
          predcure <- product
          for (kk in 1:num.batches){
            predcure[bi==kk] <- product[bi==kk]*vec.d[bi==kk]/sum(product[bi==kk])
            predcure[bi==kk] <- predcure[bi==kk] * exp(-predcure[bi==kk])
          }
          
          fitcox <- cv.fit[[i]]$fit[[j]][[k]]$fitcox
          predsurvexp <- predict(fitcox,newx=x[foldid==i,],s=min(fitcox$lambda),type="response")
          haz <- cv.fit[[i]]$fit[[j]][[k]]$haz
          cumhaz <- cv.fit[[i]]$fit[[j]][[k]]$cumhaz
          
          predcens <- survfit(Surv(ti,1-di)~1)$surv
          # tcens <- survfit(Surv(ti,1-di)~1)$time
          
          tbrier <- sort(ti)
          brier <- rep(NA,length(tbrier))
          
          for (l in 1:length(tbrier)){
            
            predsurv <- ipw <- rep(NA,length(ti)) # Conditional survival prob for all patients at tbrier[l]
            
            for (ii in 1:num.batches){
              
              bidx <- which(bi == ii)
              
              if (tbrier[l] < min(tj[[ii]])){
                
                predsurv[bidx] <- 1
                ipw[bidx] <- 1
                
              }else if (tbrier[l] >= min(tj[[ii]])){
                
                ids <- which(tj[[ii]] == max(tj[[ii]][tj[[ii]] <= tbrier[l]]))
                predsurv[bidx] <- exp(-cumhaz[[ii]][ids]*predsurvexp[bidx])
                ipw[bidx] <- as.numeric(ti[bidx] > tbrier[l])*predcens[l] + as.numeric(ti[bidx] <= tbrier[l])*predcens[bidx]
                
              }
              
            }
            
            preds <- 1 - predcure + predcure*predsurv
            brier[l] <- mean(I(ti > tbrier[l])*(1 - preds)^2/(ipw+0.001) + I(di==1 & ti <= tbrier[l])*(0 - preds)^2/(ipw+0.001))
            
          }
          
          # plot(tbrier,brier,type="S",ylim=c(0,1))
          
          cv_brier[j,k,i] <- trapz(tbrier,brier)
          
        }
        
      }
      
    }
    
  }
  
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
              foldid = foldid,
              selected = list(min=selectedmin,
                              `1se`=selected1se)
  )
  )
  
}