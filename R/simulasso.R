#' Simulate Data 
#'
#' @param n Number of patients to simulate
#' @param p Number of parameters to simulate
#' @param id_cure Which parameters to include in cure simulation
#' @param id_cox Which parameters to include in cox simulation
#' @param coefs_cure Pre-specified coefficients of cure
#' @param coefs_cox  Pre-specified Coefficients of cox
#' @param n_true_p Number of "true" non-zero coefficients
#'
#' @return
#' - `t`- simulated survival times
#' - `d`- simulated event indicator
#' - `x` simulated data n x p matrix or data.frame
#' - `id_cure` - IDs of selected cure parameters
#' - `id_cure` - IDs of selected cox parameters
#' - `coefs_cure` - Coefficients of cure parameters
#' - `coefs_cox` - Coefficients of cox parameters
#' 
#' @export
#'
#' @examples
#' simulasso(n = 10, p = 3)
simulasso <- function(n = 100,
                      p = 200,
                      id_cure = NULL,
                      id_cox = NULL,
                      coefs_cure = NULL,
                      coefs_cox = NULL, 
                      n_true_p = 10) {
  require(mvtnorm)

  t <- rep(NA, length = n)
  d <- rep(NA, length = n)

  # number of predictors is size 
  if (is.null(id_cure)) id_cure <- sample(p, size = n_true_p)
  if (is.null(id_cox)) id_cox <- sample(p, size = n_true_p)

  if (is.null(coefs_cure)) coefs_cure <- runif(length(id_cure) + 1, min = 2, max = 3)
  if (is.null(coefs_cox)) coefs_cox <- runif(length(id_cox), min = 2, max = 3)

  x <- mvtnorm::rmvnorm(n, mean = rep(0, p))

  # add intercept and calculate linear predictor and probability of cure using the logistic function:
  
  puncure <- exp(cbind(1, x[, id_cure]) %*% coefs_cure) / (1 + exp(cbind(1, x[, id_cure]) %*% coefs_cure))
  
  # make into event index
  uncure <- rbinom(n, 1, prob = puncure)

  # linear predictor survival
  lp_cox <- exp(x[, id_cox] %*% coefs_cox)

  ### Should censoring time for cured subjects be generated from a same distribution as those uncured?
  
  for(i in 1:n){
    U <- runif(1,min=0,max=1)
    t0 <- log(1-(log(1-U))/(0.1*lambda[i]))
    c0 <- min(rexp(1,rate=1/30),runif(1,min=20,max=40))
    t[i] <- ifelse(uncure[i]==1, min(t0,c0), c0)
    d[i] <- ifelse(uncure[i]==1, as.numeric(I(t0 <= c0)), 0)

  }

  dat <- list(t = t, d = d, x = x,
              uncure = uncure,
              id_cure = id_cure,
              id_cox = id_cox, coefs_cure = coefs_cure, 
              coefs_cox = coefs_cox)
}
