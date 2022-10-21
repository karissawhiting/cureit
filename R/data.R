#' Results from a simulated study of two chemotherapy agents
#'
#' A dataset containing the baseline characteristics of 200 patients
#' who received Drug A or Drug B.  Dataset also contains the outcome of
#' tumor response to the treatment.
#'
#' @format A data frame with 200 rows--one row per patient
#' \describe{
#'     \item{trt}{Chemotherapy Treatment}
#'     \item{age}{Age}
#'     \item{marker}{Marker Level (ng/mL)}
#'     \item{stage}{T Stage}
#'     \item{grade}{Grade}
#'     \item{response}{Tumor Response}
#'     \item{death}{Patient Died}
#'     \item{ttdeath}{Months to Death/Censor}
#' }
"trial"

#' Survival after malignant melanoma
#'
#' Data frame has 205 rows and 7 columns.
#' It contains clinical variables and survival outcome data relating to the survival of patients after an operation for
#' malignant melanoma. It was collected at Odense University Hospital by K.T. Drzewiecki and originally sourced from the {ISwR} Package
#'
#'
#' @format This data frame contains the following columns:
#'
#' \describe{
#'       \item{id}{a numeric vector, patient code.}
#'       \item{status}{a numeric vector code, survival status; 1: dead from melanoma, 0: alive or dead from other cause.}
#'       \item{days}{a numeric vector, observation time in days.}
#'       \item{years}{a numeric vector, observation time in years.}
#'       \item{ulceration}{a numeric vector code, ulceration; 1: present, 2: absent.}
#'       \item{tumor_thickness}{a numeric vector, tumor thickness (1/100 mm).}
#'       \item{sex}{a numeric vector code; 1: female, 2: male.}
#'     }
#' @source P.K. Andersen, \enc{Ø}{O}. Borgan, R.D. Gill, and N. Keiding (1991),
#'   \emph{Statistical Models Based on
#'     Counting Processes}, Appendix 1, Springer-Verlag.
"melanoma"
