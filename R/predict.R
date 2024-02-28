#' #' Predicted survival probability by cure model
#' #'
#' #' @param object A cureit object.
#' #' @param times Numeric vector of times to obtain survival probability estimates at
#' #' @param probs Numeric vector of quantiles to obtain estimates at
#' #' @param newdata A `base::data.frame()` or `tibble::tibble()` containing all
#' #' the original predictors used to create object. Defaults to `NULL`.
#' #' @param type Output format of predicted values: "lp" (linear predictor) or "prob" (predicted probabilities).
#' #' @param brier Boolean indicator of calculating the Brier scores at specified `times`.
#' #' @param cox Boolean indicator of fitting the Cox model for the training data and calculating the Brier scores at specified `times` for `newdata`.
#' #' @inheritParams cureit
#' #' @return named list of prediction estimates
#' #' @family cureit() functions
#' #' @export
#' #' @examples
#' #' p <- cureit(surv_formula = Surv(ttdeath, death) ~ age, 
#' #'    cure_formula = ~ age,
#' #'    data = trial) 
#' #'    
#' #' pred <- predict(p, times = seq(5,24,0.5),
#' #'  newdata = trial[complete.cases(trial), ], brier=TRUE,cox=TRUE)
#' #' 
#' #' plot(seq(5,24,0.5),pred$brier,type="S",pch=1)
#' #' lines(seq(5,24,0.5),pred$brier_cox,type="S",col="red",pch=3)
#' #' legend("topright",c("Cure model","Cox model"),
#' #' col=c("black","red"),lty=1,pch=c(1,3))
#' #' 
#' #'   
#' predict.cureit <- function(object, times = NULL, probs = NULL,
#'                            newdata = NULL, type = c("prob", "lp"),
#'                            brier = FALSE, cox = FALSE, ...) {
#'   # checking inputs ------------------------------------------------------------
#'   
#'   type <- match.arg(type)
#'   
#'   if (is.null(times) + is.null(probs) != 1L) {
#'     stop("Must specify one and only one of `times=` and `probs=`.", call. = FALSE)
#'   }
#'   if (!is.null(times) && any(times < 0)) {
#'     stop("`times=` must be non-negative.", call. = FALSE)
#'   }
#'   if (!is.null(probs) && !all(dplyr::between(probs, 0, 1))) {
#'     stop("`probs=` must be between 0 and 1.", call. = FALSE)
#'   }
#'   if (!(type %in% c("lp","prob"))) {
#'     stop("`type=` must be one out of 'lp' or 'prob'.", call. = FALSE)
#'   }
#'   # if (type == "lp" & brier){
#'   #   stop("`type=` must be 'prob' to calculate the Brier scores.")
#'   # }
#'   # if (brier & is.null(times)){
#'   #   stop("Must specify one or more times to calculate the Brier scores.")
#'   # }
#'   # if (!brier & cox){
#'   #   stop("Must let `brier=TRUE` to perform Cox regression and corresponding Brier scores.")
#'   # }
#'   
#'   # getting predictions on the original model fit ------------------------------
#'   processed <- cureit_mold(object$surv_formula, object$cure_formula, newdata %||% object$data,
#'                            surv_blueprint = object$surv_blueprint, cure_blueprint = object$cure_blueprint)
#'   
#'   newX=processed$surv_processed$predictors
#'   newZ=processed$cure_processed$predictors
#'   if (is.vector(newZ)) 
#'     newZ = as.matrix(newZ)
#'   newZ = cbind(1, newZ)
#'   if (is.vector(newX)) 
#'     newX = as.matrix(newX)
#'   
#'   s0 = as.matrix(object$smcure$s[order(object$smcure$Time)], ncol = 1)
#'   n = nrow(s0)
#'   uncureprob = exp(object$smcure$b %*% t(newZ))/(1 + exp(object$smcure$b %*% t(newZ)))
#'   
#'   scure = array(0, dim = c(n, nrow(newX)))
#'   t = array(0, dim = c(n, nrow(newX)))
#'   spop = array(0, dim = c(n, nrow(newX)))
#'   
#'   #### Only support PH model for now
#'   ebetaX = exp(object$smcure$beta %*% t(newX))
#'   for (i in 1:nrow(newZ)) {
#'     scure[, i] = s0^ebetaX[i]
#'   }
#'   
#'   for (i in 1:n) {
#'     for (j in 1:nrow(newX)) {
#'       spop[i, j] = uncureprob[j] * scure[i, j] + (1 - 
#'                                                     uncureprob[j])
#'     }
#'   }
#'   
#'   spop_prd = cbind(object$smcure$Time, spop)
#'   scure_prd = cbind(object$smcure$Time, scure)
#'   cure_prd = 1-uncureprob
#'   
#'   # Type = "Prob ---------
#'   if (type == "prob"){
#'     
#'     if (!is.null(times)) {
#'       return(list(cured = cure_prd,
#'                     surv_uncured = probs_at_times(scure_prd,times),
#'                     surv_marginal = probs_at_times(spop_prd,times)))
#'     }
#'     
#'     if (!is.null(probs)) {
#'       return(list(cured = cure_prd,
#'                   surv_uncured = times_at_probs(scure_prd,probs),
#'                   surv_marginal = times_at_probs(spop_prd,probs)))  
#'     }
#'   }
#'   
#'   # Type = "Prob ---------
#'   if (type == "lp"){
#'     
#'     return(list(lp_cure_model = object$smcure$b %*% t(newZ),
#'                 lp_surv_model = object$smcure$beta %*% t(newX)))
#'     
#'   }
#'   
#' }
#' 
#' times_at_probs <- function(matrix_pred, probs) {
#'   matrix_zero <- matrix(c(0, 1), ncol = 2)
#'   lst_time_risk <-
#'     purrr::map(
#'       probs,
#'       ~ purrr::map_dbl(
#'         seq_len(ncol(matrix_pred) - 1L),
#'         function(i) {
#'           # adding 1 to the matrix, keeping only the time column, and col of interest
#'           m <- rbind(matrix_zero, matrix_pred[, c(1L, i + 1L)])
#'           
#'           # return NA if the quantiles are all missing OR prob is larger than observed
#'           if (isTRUE(all(is.na(m[-1, 2])) || .x > max(m[, 2]))) {
#'             return(NA)
#'           }
#'           if (isTRUE(.x %in% m[, 2])) {
#'             return(m[m[, 2] %in% .x, 1])
#'           }
#'           return(m[which.min(m[, 2] >= .x), 1])
#'         }
#'       )
#'     ) %>%
#'     stats::setNames(paste0("prob ", probs * 100, "%"))
#'   
#'   # returning results ----------------------------------------------------------
#'   lst_time_risk
#' }
#' 
#' 
#' 
#' probs_at_times <- function(matrix_pred, times) {
#'   # defining times for predictions ---------------------------------------------
#'   # CHECK- here should this be warning?
#'   all_times <- union(0, matrix_pred[, 1]) %>% sort()
#'   if (isTRUE(max(times) > max(all_times))) {
#'     stringr::str_glue("`times=` cannot be larger than {max(all_times)}") %>%
#'       stop(call. = FALSE)
#'   }
#'   times_obs <-
#'     purrr::map_dbl(
#'       times,
#'       function(.x) {
#'         if (isTRUE(.x %in% all_times)) {
#'           return(all_times[all_times %in% .x])
#'         }
#'         all_times[which.min(all_times <= .x) - 1L]
#'       }
#'     )
#'   
#'   # named list of the risks, the names are the times,
#'   # the values are the estimates of risk at the covar levels
#'   lst_risk_time <-
#'     purrr::map(seq_len(length(all_times) - 1L), ~ matrix_pred[.x, -1]) %>%
#'     stats::setNames(all_times[-1]) %>%
#'     dplyr::bind_cols() %>%
#'     mutate(`0` = 0, .before = 1) %>%
#'     as.list() %>%
#'     stats::setNames(all_times)
#'   
#'   # extracting risks at specified time -----------------------------------------
#'   lst_risk_time[as.character(times_obs)] %>%
#'     stats::setNames(paste("time", times))
#' }
