#' 
#' 
#' 
#' mod <-
#'   coxph(Surv(rfs_years, local_recurrence) ~
#'           primary_tumor_location +
#'           primary_tumor_size_cat +
#'           mutations_grouped,
#'         data = msk2, x =  TRUE)
#' 
#' risk <- predict(mod, type = "risk")   # risk from Cox PH model
#' 
#' 
#' #' K-index estimation
#' #'
#' #' @export
#' #' @param risk the risk score for survival part
#' #' @param pi the uncured probability
#' #' @param model the type of survival model ("PH", "PO","Normal")
#' 
#' 
#' k_index <- function(object, times,...) {
#'   
#'   # risk score for survival
#'   est.risk <- X %*% beta
#'   
#'   # uncured probability 
#'   est.odds <- Z %*% b
#'   est.pi   <- logit.inv(est.odds)
#'   N <- sum(pi)
#'   comp.risk.cox <- outer(risk, risk, .fun, model = model)
#'   diag(comp.risk.cox) <- 0
#'   k <- sum(comp.risk.cox * outer(pi, pi, "*")  ) * 2  / ( N^2 - sum(pi^2)) # Based on definition
#'   k
#'   
#'   if(toupper(model) == "PH" & baseline == T){ 
#'     est.surv <- surv ^ exp( est.risk) 
#'   }else{
#'       est.surv <- surv}
#'   
#'   est.w <- w.cure(est.pi, delta, est.surv)
#'   
#'   #   risk <- as.numeric(risk)
#'   #   pi   <- as.numeric(pi)
#'   #   N <- sum(pi)
#'   #   comp.risk.cox <- outer( risk, risk, .fun, model = model)
#'   #   diag(comp.risk.cox) <- 0
#'   #   k <- sum(comp.risk.cox *   outer(pi, pi, "*")  ) * 2  / ( N^2 - sum(pi^2)) # Based on definition
#'   #   k
#'   
#' }
#' 
#' 
#'   # Data checks -------
#'   if (object$nboot <= 0) {
#'     stop("`nboot=` cannot be less than 1. Please rerun `cureit` with `nboot` > 0", call. = FALSE)
#'   }
#'   if (!is.null(times) && any(times < 0)) {
#'     stop("`times=` must be non-negative.", call. = FALSE)
#'   }
#'   
#'   original <- predict(
#'     object = object,
#'     newdata = object$data,
#'     method = "prob",
#'     times = times, 
#'     brier = TRUE, cox = TRUE
#'   )
#'   
#' }
#' 
#' # test --------------------------------------------------------------------
#' logit.inv <- function(x) 1/ (1 + exp(-x))
#' 
#' 
#' trial <- na.omit(trial)
#' p <- cureit(surv_formula = Surv(ttdeath, death) ~ age,
#'    cure_formula = ~ age,
#'    data = trial)
#' 
#' fit <- p$smcure
#' 
#' fit$beta
#' 
#' # X 
#' mf <- model.frame(p$surv_formula, p$data)
#' n <- dim(data)[1]
#' 
#' cvars <- all.vars(p$cure_formula)
#' Z <- as.matrix(cbind(rep(1,n),data[,cvars]))
#' colnames(Z) <- c("(Intercept)",cvars)
#' 
#' offsetvar <- NULL
#' Y <- model.extract(mf,"response")
#' X <- model.matrix(attr(mf,"terms"), mf)
#' X <- X[, 2]
#' 
#' Time <- Y[,1]
#' Status <- Y[,2]
#' 
#' 
#' Z <- Z
#' time <- Time
#' delta <- Status
#' beta <- fit$beta
#' b <- fit$b
#' surv <- fit$s
#' 
#' 
#' # risk score for survival
#' est.risk <- X * beta
#' 
#' # uncured probability 
#' est.odds <- Z * b
#' est.pi   <- logit.inv(est.odds)
#' 
#' N <- sum(est.pi)
#' comp.risk.cox <- outer(est.risk, est.risk, .fun, model = "PH")
#' diag(comp.risk.cox) <- 0
#' k <- sum(comp.risk.cox * outer(est.pi, est.pi, "*")  ) * 2  / ( N^2 - sum(est.pi^2)) # Based on definition
#' k
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' eva_cure(time,delta,X,beta,Z,b,surv,model)
#' 
#' .fun <- function(x1,x2, model){  # Kernel Function
#'   diff <- x1 - x2
#'   
#'   if( toupper(model) == "PH") weight <- 1/(1+exp(diff))
#'   if( toupper(model) == "PO"){ weight <- (1 - exp(diff) + diff * exp(diff)) / (1 - exp(diff))^2
#'   weight <- ifelse(is.nan(weight), 0.5, weight)
#'   }
#'   if( toupper(model) == "NORMAL") weight <- 1 - pnorm(diff/ sqrt(2) )
#'   (diff < 0) * weight + 0.5 * (diff == 0)
#' }
#' 
#' k.ind <- function(risk, pi, model){
#'   risk <- as.numeric(risk)
#'   pi   <- as.numeric(pi)
#'   N <- sum(pi)
#'   comp.risk.cox <- outer( risk, risk, .fun, model = model)
#'   diag(comp.risk.cox) <- 0
#'   k <- sum(comp.risk.cox *   outer(pi, pi, "*")  ) * 2  / ( N^2 - sum(pi^2)) # Based on definition
#'   k
#' }
#' 
#' 
