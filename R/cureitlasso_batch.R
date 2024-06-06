cureitlasso_batch <- function(t,
                              d,
                              x,
                              b, # batch labels
                              minmu.ratio=0.01, # minimum penalty value in the logistic model
                              minlambda.ratio=0.01, # minimum penalty value in the cox model
                              adaptive=FALSE,
                              length.grid=10,
                              mus=NULL,
                              lambdas=NULL,
                              tol=1e-2,
                              maxiter=100,
                              progress=FALSE
){
  
  require(glmnet)
  require(survival)
  require(clogitL1)
  
  tau <- matrix(0.5,nrow=length(t),ncol=2) # First column: cured; second column: uncured
  tau[d==1,1] = 0
  tau[d==1,2] = 1
  lab_class = c(rep(0,length(t)),rep(1,length(t))) # Pseudo class labels (0=cured, 1=uncured)
  
  tj <- sort(t[d==1])
  
  # CV: Be careful for fold split! Two replicates from the same subject should always be assigned to the same fold
  
  penalty.factor.cure=rep(1,ncol(x))
  penalty.factor.cox=rep(1,ncol(x))
  
  # if (!adaptive){
  #   
  # }else{
  #   
  #   ridgecure <- cv.glmnet(x=do.call("rbind", rep(list(x), 2)),
  #                          y=lab_class,
  #                          family="binomial",
  #                          weights=as.vector(tau),
  #                          lambda.min.ratio = minmu.ratio,
  #                          alpha=0)
  #   
  #   penalty.factor.cure <- 1/abs(coef(ridgecure, s=ridgecure$lambda.min)[-1,1])
  #   
  #   ridgecox <- cv.glmnet(x=x,
  #                         y=Surv(t,d),
  #                         family="cox",
  #                         weights=tau[,2],
  #                         lambda.min.ratio = minlambda.ratio,
  #                         alpha=0)
  #   
  #   penalty.factor.cox <- 1/abs(coef(ridgecox, s=ridgecox$lambda.min)[,1])
  #   
  # }
  
  num.batches <- length(unique(b))
  # num.d <- rep(NA,num.batches)
  vec.d <- rep(NA,nrow(x))
  
  for (i in 1:num.batches){
    vec.d[b==i] <- sum(d[b == i])
  }
  
  fitcure0 <- clogitL1(x, d, b,
                       minLambdaRatio = minmu.ratio)
  
  product <- exp(x %*% t(fitcure0$beta[-1,])) # rows: samples; columns: penalty parameters
  pred <- product
  for (i in 1:num.batches){
    pred[b==i,] <- t(apply(product[b==i,]*vec.d[b==i], 1, function(x) x/colSums(product[b==i,])))
    pred[b==i,] <- pred[b==i,] * exp(-pred[b==i,])
  }
  
  fitcox0 <- glmnet(x=x,
                    y=stratifySurv(Surv(t,d),b),
                    family="cox",
                    weights=tau[,2],
                    lambda.min.ratio = minlambda.ratio,
                    penalty.factor = penalty.factor.cox)
  
  if (is.null(mus) | is.null(lambdas)){
    idx_mus <- seq(2,length(fitcure0$lambda),length.out=length.grid)
    mus <- fitcure0$lambda[idx_mus]
    idx_lambdas <- seq(1,length(fitcox0$lambda),length.out=length.grid)
    lambdas <- fitcox0$lambda[idx_lambdas]
  }
  
  musmax <- fitcure0$lambda[2]
  lambdasmax <- fitcox0$lambda[1]
  
  predcure <- pred #Prob of uncured
  predsurvexp <- predict(fitcox0,newx=x,s=lambdas,type="response") # exp(betaX)
  
  tj <- list()
  haz <- list()
  cumhaz <- list()
  
  for (j in 1:num.batches){
    
    tj[[j]] <- sort(t[d==1 & b==j])
    haz[[j]] <- matrix(0,nrow=length(tj[[j]]),ncol=length(lambdas))
    for (i in 1:length(tj[[j]])){
      haz[[j]][i,] <- 1/colSums(tau[t >= tj[[j]][i] & b==j,2] * predsurvexp[t >= tj[[j]][i] & b==j,])
    }
    
    cumhaz[[j]] <- apply(haz[[j]],2,cumsum)
    
  }
  
  predsurv <- matrix(0,nrow=length(t),ncol=length(lambdas))
  for (i in 1:length(t)){
    
    for (j in 1:num.batches){
      
      if (b[i] == j){
        
        if (t[i] %in% tj[[j]]){
          predsurv[i,] <- haz[[j]][tj[[j]]==t[i],]*predsurvexp[i,]*exp(-cumhaz[[j]][tj[[j]]==t[i],]*predsurvexp[i,])
        }else if (t[i] > min(tj[[j]])){
          ids <- which(tj[[j]] == max(tj[[j]][tj[[j]] < t[i]]))
          predsurv[i,] <- exp(-cumhaz[[j]][ids,]*predsurvexp[i,])
        }
        
      }
      
    }
    
  }
  
  ### Given specific i and j: then do iterations
  
  fit <- list()
  # num_alpha <- matrix(NA,nrow=length(mus),ncol=length(lambdas))
  # num_beta <- matrix(NA,nrow=length(mus),ncol=length(lambdas))
  
  for (i in 1:length(mus)){
    
    fit[[i]] <- list()
    
    for (j in 1:length(lambdas)){
      
      tau <- d + (1-d) * (predcure[,i] * predsurv[,j])/( (1 - predcure[,i]) + predcure[,i] * predsurv[,j] )
      loglik <- sum((1-d)*log(1-predcure[,i] + predcure[,i] * predsurv[,j] + 0.001) + d*log(predcure[,i] * predsurv[,j] + 0.001))
      
      if (mus[i] > musmax){
        alpha <- fitcure0$beta[max(which(fitcure0$lambda >= musmax)),]
      }else{
        alpha <- fitcure0$beta[max(which(fitcure0$lambda >= mus[i])),]
      }
      
      if (lambdas[j] > lambdasmax){
        beta <- fitcox0$beta[,max(which(fitcox0$lambda >= lambdasmax))]
      }else{
        beta <- fitcox0$beta[,max(which(fitcox0$lambda >= lambdas[j]))]
      }
      
      for (iter in 1:maxiter){
        
        tau0 <- tau
        alpha0 <- alpha
        beta0 <- beta
        loglik0 <- loglik
        
        # if (!adaptive){
        #   
        #   penalty.factor.cure=rep(1,ncol(x))
        #   penalty.factor.cox=rep(1,ncol(x))
        #   
        # }else{
        #   
        #   ridgecure <- cv.glmnet(x=do.call("rbind", rep(list(x), 2)),
        #                          y=lab_class,
        #                          family="binomial",
        #                          weights=c(1-tau0,tau0),
        #                          lambda.min.ratio = minmu.ratio,
        #                          alpha=0)
        #   
        #   penalty.factor.cure <- 1/abs(coef(ridgecure, s=ridgecure$lambda.min)[-1,1])
        #   
        #   ridgecox <- cv.glmnet(x=x,
        #                         y=Surv(t,d),
        #                         family="cox",
        #                         weights=tau0,
        #                         lambda.min.ratio = minlambda.ratio,
        #                         alpha=0)
        #   
        #   penalty.factor.cox <- 1/abs(coef(ridgecox, s=ridgecox$lambda.min)[,1])
        #   
        # }
        
        vec.d <- rep(NA,nrow(x))
        
        for (k in 1:num.batches){
          vec.d[b==k] <- sum(tau0[b==k] >= 0.5)
        }
        
        fitcure0 <- clogitL1(x, 
                             as.numeric(tau0 >= 0.5), 
                             b,
                             numLambda = 100)
        
        fitcure <- clogitL1(x,
                            as.numeric(tau0 >= 0.5), 
                            b,
                            minLambdaRatio = mus[i]/fitcure0$lambda[1])
        
        product <- exp(x %*% fitcure$beta[length(fitcure$lambda),]) # rows: samples; columns: penalty parameters
        pred <- product
        for (k in 1:num.batches){
          pred[b==k] <- product[b==k]*vec.d[b==k]/sum(product[b==k])
          pred[b==k] <- pred[b==k] * exp(-pred[b==k])
        }
        
        fitcox <- glmnet(x=x,
                         y=stratifySurv(Surv(t,d),b),
                         family="cox",
                         weights=tau0,
                         lambda=10^(seq(log10(lambdasmax*10),log10(lambdas[j]),length.out=100)),
                         penalty.factor = penalty.factor.cox)
        
        predcure_iter <- pred #Prob of uncured
        predsurvexp <- predict(fitcox,newx=x,s=fitcox$lambda[length(fitcox$lambda)],type="response") # exp(betaX)
        
        tj <- list()
        haz <- list()
        cumhaz <- list()
        
        for (k in 1:num.batches){
          
          tj[[k]] <- sort(t[d==1 & b==k])
          haz[[k]] <- rep(0,nrow=length(tj[[k]]))
          for (h in 1:length(tj[[k]])){
            haz[[k]][h] <- 1/sum(tau0[t >= tj[[k]][h] & b==k] * predsurvexp[t >= tj[[k]][h] & b==k,])
          }
          
          cumhaz[[k]] <- cumsum(haz[[k]])
          
        }
        
        predsurv_iter <- rep(0,length(t))
        for (ii in 1:length(t)){
          
          for (jj in 1:num.batches){
            
            if (b[ii] == jj){
              
              if (t[ii] %in% tj[[jj]]){
                predsurv_iter[ii] <- haz[[jj]][tj[[jj]]==t[ii]]*predsurvexp[ii,]*exp(-cumhaz[[jj]][tj[[jj]]==t[ii]]*predsurvexp[ii,])
              }else if (t[ii] > min(tj[[jj]])){
                ids <- which(tj[[jj]] == max(tj[[jj]][tj[[jj]] < t[ii]]))
                predsurv_iter[ii] <- exp(-cumhaz[[jj]][ids]*predsurvexp[ii,])
              }
              
            }
            
          }
          
        }
        
        tau <- d + (1-d) * (predcure_iter * predsurv_iter)/( (1 - predcure_iter) + predcure_iter * predsurv_iter )
        alpha <- fitcure$beta[length(fitcure$lambda),]
        beta <- fitcox$beta[,length(fitcox$lambda)]
        
        difftau <- mean((tau-tau0)^2)
        diffpar <- max(abs(alpha-alpha0),abs(beta-beta0))
        loglik <- sum((1-d)*log(1-predcure_iter + predcure_iter * predsurv_iter + 0.001) + d*log(predcure_iter * predsurv_iter + 0.001))
        
        if(iter == 1){
          l0 = loglik
        }else if(iter == 2){
          diffl = loglik - l0
          l0 = loglik
        }else if(iter == 3){
          a = (loglik - l0)/(diffl)
          diffl = loglik - l0
          lA = l0 + diffl/(1-a)
          l0 = loglik
        }else if(iter >= 4){
          a = (loglik - l0)/(diffl)
          diffl = loglik - l0
          difflA = l0 + diffl/(1-a) - lA
          lA = l0 + diffl/(1-a)
          l0 = loglik
          
          if (abs(difflA) < tol){
            
            if (progress){
              print(i)
              print(j)
              print(iter)
              print(paste0("difflA:",difflA))
              print(paste0("diffpar:",diffpar))
              print(paste0("difftau:",difftau))
              print(paste0("difflik:",loglik-loglik0))
            }
            
            # if (adaptive){
            #   
            #   ridgecure <- cv.glmnet(x=do.call("rbind", rep(list(x), 2)),
            #                          y=lab_class,
            #                          family="binomial",
            #                          weights=c(1-tau,tau),
            #                          lambda.min.ratio = minmu.ratio,
            #                          alpha=0)
            #   
            #   penalty.factor.cure <- 1/abs(coef(ridgecure, s=ridgecure$lambda.min)[-1,1])
            #   
            #   ridgecox <- cv.glmnet(x=x,
            #                         y=Surv(t,d),
            #                         family="cox",
            #                         weights=tau,
            #                         lambda.min.ratio = minlambda.ratio,
            #                         alpha=0)
            #   
            #   penalty.factor.cox <- 1/abs(coef(ridgecox, s=ridgecox$lambda.min)[,1])
            #   
            #   fitcure <- glmnet(x=do.call("rbind", rep(list(x), 2)),
            #                     y=lab_class,
            #                     family="binomial",
            #                     weights=c(1-tau,tau),
            #                     lambda=10^(seq(musmax*10,log10(mus[i]),length.out=100)),
            #                     penalty.factor = penalty.factor.cure)
            #   
            #   fitcox <- glmnet(x=x,
            #                    y=Surv(t,d),
            #                    family="cox",
            #                    weights=tau,
            #                    lambda=10^(seq(lambdasmax*10,log10(lambdas[i]),length.out=100)),
            #                    penalty.factor = penalty.factor.cox)
            #   
            #   predcure_iter <- predict(fitcure,newx=x,s=min(fitcure$lambda),type="response") #Prob of uncured
            #   predsurvexp_iter <- predict(fitcox,newx=x,s=min(fitcox$lambda),type="response") # exp(betaX)
            #   
            #   haz <- rep(0,nrow=length(tj))
            #   for (k in 1:length(tj)){
            #     haz[k] <- 1/sum(tau0[t >= tj[k]] * predsurvexp_iter[t >= tj[k]])
            #   }
            #   cumhaz <- cumsum(haz)
            #   
            #   predsurv_iter <- rep(0,nrow=length(t))
            #   for (k in 1:length(t)){
            #     if (t[k] %in% tj){
            #       predsurv_iter[k] <- haz[tj==t[k]]*predsurvexp_iter[k]*exp(-cumhaz[tj==t[k]]*predsurvexp_iter[k])
            #     }else if (t[k] > min(tj)){
            #       ids <- which(tj == max(tj[tj < t[k]]))
            #       predsurv_iter[k] <- exp(-cumhaz[ids]*predsurvexp_iter[k])
            #     }
            #   }
            #   
            #   tau <- d + (1-d) * (predcure_iter * predsurv_iter)/( (1 - predcure_iter) + predcure_iter * predsurv_iter )
            #   alpha <- fitcure$beta[,length(fitcure$lambda)]
            #   beta <- fitcox$beta[,length(fitcox$lambda)]
            #   
            # }
            
            break
          } 
        }
        
      }
      
      fit[[i]][[j]] <- list(t=t,
                            tj=tj,
                            d=d,
                            x=x,
                            b=b,
                            # alpha0=fitcure$a0[length(fitcure$lambda)],
                            alpha=alpha,
                            beta=beta,
                            haz=haz,
                            cumhaz=cumhaz,
                            tau=as.vector(tau),
                            predsurv=as.vector(predsurv_iter),
                            predcure=as.vector(predcure_iter),
                            fitcure=fitcure,
                            fitcox=fitcox)
      
      # num_alpha[i,j] <- sum(alpha!=0)
      # num_beta[i,j] <- sum(beta!=0)
      
    }
    
  }
  
  return(list(fit=fit,
              # num_alpha=num_alpha,
              # num_beta=num_beta,
              mus=mus,
              lambdas=lambdas))
  
}