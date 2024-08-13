#' gtsummary methods
#'
#' @description
#' This function is an S3 methods for working with [`cureit::cureit()`] model
#' results.
#'
#' - `tbl_regression.cureit()`: This function sets the cureit tidier for `cureit()` models.
#' @param x (`cureit`)\cr
#'   `cureit::cureit()` regression object
#' @param tidy_fun (`function`)\cr
#'    Tidier function for the model. Default is `cureit::tidy()`.
#' @param type not used
#' @inheritParams gtsummary::tbl_regression
#'
#' @return gtsummary table or data frame of results
#' @export
#'
#' @examples
#' cureit(surv_formula = Surv(days, status) ~ ulceration,
#'   cure_formula = ~ ulceration, data = melanoma) |>
#'   gtsummary::tbl_regression() |>
#'   gtsummary::as_gt()
#'   
tbl_regression.cureit <- function(x, tidy_fun = cureit::tidy, ...) {
  suppressMessages(
    asNamespace("gtsummary")[["tbl_regression.default"]](x = x, tidy_fun = tidy_fun, ...))
}
