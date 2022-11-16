#' Inference for Brier score for cureit objects
#' @rdname Brier_inference_cureit
#' @param object A cureit object.
#' @param newdata A `base::data.frame()` or `tibble::tibble()` containing all
#' the original predictors used to create x. Defaults to `NULL`.
#' @param nboot Number of bootstrap resamples generated.
#' @param times Numeric vector of times to obtain survival probability estimates at
#' @param ... Additional arguments passed to other methods.
#' @return a tibble
#' @family cureit() functions
#' @export
#' @examples
#' trial <- na.omit(trial)
#' p <- cureit(
#'   surv_formula = Surv(ttdeath, death) ~ age,
#'   cure_formula = ~age,
#'   data = trial
#' )
#'
#' times <- seq(5, 24, 0.5)
#' bootbrier <- Brier_inference(p, times = times)
#' plot(times, bootbrier$brier, type = "s", xlab = "Time", ylab = "Brier score")
#' lines(times, bootbrier$brier_2.5, type = "s", lty = 2)
#' lines(times, bootbrier$brier_97.5, type = "s", lty = 2)
#' lines(times, bootbrier$brier_cox, type = "s", col = "red")
#' lines(times, bootbrier$brier_cox_2.5, type = "s", lty = 2, col = "red")
#' lines(times, bootbrier$brier_cox_97.5, type = "s", lty = 2, col = "red")
#' legend("topleft", c("Cure model", "Cox model"), col = c("black", "red"), lty = 1)
#'
Brier_inference <- function(object,
                            newdata = NULL,
                            nboot = 100,
                            times,
                            ...) {

  # Data checks -------
  if (nboot <= 0) {
    stop("`nboot=` cannot be less than 1", call. = FALSE)
  }
  if (!is.null(times) && any(times < 0)) {
    stop("`times=` must be non-negative.", call. = FALSE)
  }

  Brier_original <- predict(
    object = object,
    newdata = newdata, method = "prob", times = times, brier = TRUE, cox = TRUE
  )

  newdata <- newdata %||% object$data
  processed <- cureit_mold(object$surv_formula, object$cure_formula, newdata %||% object$data)
  s.outcomes <- as.matrix(processed$surv_processed$outcomes[, 1, drop = TRUE])
  status <- s.outcomes[, "status"]

  data1 <- subset(newdata, status == 1)
  data0 <- subset(newdata, status == 0)
  n1 <- nrow(data1)
  n0 <- nrow(data0)

  boot_brier <- boot_brier_cox <- matrix(NA, nrow = nboot, ncol = length(times))

  # Start bootstrapping

  i <- 1

  while (i <= nboot) {

    # print(i)

    id1 <- sample(1:n1, n1, replace = TRUE)
    id0 <- sample(1:n0, n0, replace = TRUE)
    bootdata <- rbind(data1[id1, ], data0[id0, ])

    bootfit <- cureit(
      surv_formula = object$surv_formula,
      cure_formula = object$cure_formula,
      data = bootdata,
      nboot = 0,
      eps = object$eps
    )

    pred_bootfit <- predict(bootfit, method = "prob", times = times, brier = TRUE, cox = TRUE)

    boot_brier[i, ] <- pred_bootfit$brier
    boot_brier_cox[i, ] <- pred_bootfit$brier_cox

    i <- i + 1
  }

  brier_sd <- apply(boot_brier, 2, stats::sd)
  brier_cox_sd <- apply(boot_brier_cox, 2, stats::sd)
  brier_2.5 <- apply(boot_brier, 2, function(x) stats::quantile(x, probs = 0.025))
  brier_97.5 <- apply(boot_brier, 2, function(x) stats::quantile(x, probs = 0.975))
  brier_cox_2.5 <- apply(boot_brier_cox, 2, function(x) stats::quantile(x, probs = 0.025))
  brier_cox_97.5 <- apply(boot_brier_cox, 2, function(x) stats::quantile(x, probs = 0.975))

  return(list(
    brier = Brier_original$brier,
    brier_cox = Brier_original$brier_cox,
    brier_sd = brier_sd,
    brier_cox_sd = brier_cox_sd,
    brier_2.5 = brier_2.5,
    brier_97.5 = brier_97.5,
    brier_cox_2.5 = brier_cox_2.5,
    brier_cox_97.5 = brier_cox_97.5,
    boot_brier = boot_brier,
    boot_brier_cox = boot_brier_cox
  ))
}