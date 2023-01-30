
#' K-index estimation
#'
#' @export
#' @param risk the risk score for survival part
#' @param pi the uncured probability
#' @param model the type of survival model ("PH", "PO","Normal")



# mod <- 
#   coxph(Surv(rfs_years, local_recurrence) ~ 
#           primary_tumor_location +
#           primary_tumor_size_cat +
#           mutations_grouped, 
#         data = msk2, x =  TRUE)
# 
# risk <- predict(mod, type = "risk")   # risk from Cox PH model
# 
# k_index <- function(risk, pi, model){
#   risk <- as.numeric(risk)
#   pi   <- as.numeric(pi)
#   N <- sum(pi)
#   comp.risk.cox <- outer( risk, risk, .fun, model = model)
#   diag(comp.risk.cox) <- 0
#   k <- sum(comp.risk.cox *   outer(pi, pi, "*")  ) * 2  / ( N^2 - sum(pi^2)) # Based on definition
#   k
# }