#' Logit transformation
#' @export
#' @noRd
#' @keywords internal
.kernel_fun <- function(x1,x2, model){
  diff <- x1 - x2
  weight <- 1/(1+exp(diff))
  (diff < 0) * weight + 0.5 * (diff == 0)
}

#' Logit transformation
#' @export
#' @noRd
#' @keywords internal
.logit <- function(x) log(x / (1 - x))

#' Inverse logit transformation
#' @export
#' @noRd
#' @keywords internal
.logit_inv <- function(x) 1/ (1 + exp(-x))
