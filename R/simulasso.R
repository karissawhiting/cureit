simulasso <- function(n=100,
                      p=200,
                      id_cure=NULL,
                      id_cox=NULL,
                      coefs_cure=NULL,
                      coefs_cox=NULL
                      ){
  
  require(mvtnorm)
  
  t <- rep(NA,length=n)
  d <- rep(NA,length=n)
  
  if(is.null(id_cure)) id_cure = sample(p,size=3)
  if(is.null(id_cox)) id_cox = sample(p,size=3)
  
  if(is.null(coefs_cure)) coefs_cure = runif(length(id_cure)+1,min=2,max=3)
  if(is.null(coefs_cox)) coefs_cox = runif(length(id_cox),min=2,max=3)
  
  x <- mvtnorm::rmvnorm(n,mean=rep(0,p))
  puncure <- exp(cbind(1,x[,id_cure]) %*% coefs_cure)/(1+exp(cbind(1,x[,id_cure]) %*% coefs_cure))
  uncure <- rbinom(n,1,prob=puncure)
  
  lambda <- exp(x[,id_cox] %*% coefs_cox)
  
  ### Should censoring time for cured subjects be generated from a same distribution as those uncured?
  
  for(i in 1:n){
    U <- runif(1,min=0,max=1)
    t0 <- log(1-(log(1-U))/(0.1*lambda[i]))
    c0 <- min(rexp(1,rate=0.1),runif(1,min=5,max=6))
    t[i] <- ifelse(uncure[i]==1, min(t0,c0), c0)
    d[i] <- ifelse(uncure[i]==1, as.numeric(I(t0 <= c0)), 0)
  }
  
  dat <- list(t=t,d=d,x=x,uncure=uncure,id_cure=id_cure,id_cox=id_cox,coefs_cure=coefs_cure,coefs_cox=coefs_cox)
  
}