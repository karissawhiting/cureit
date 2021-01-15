
# Clean Model Results ----------------------------------------------------------
mmcure <- function(fit_k2_cat) {
  enframe(fit_k2_cat$b) %>% 
  left_join(enframe(fit_k2_cat$b_pvalue),
            by = "name") %>%
  transmute(name, 
            value = round(exp(value.x), 2), 
            pvalue = round(value.y, 2))
}

mmsurv <- function(fit_k2_cat) {
  enframe(fit_k2_cat$beta) %>% 
  left_join(enframe(fit_k2_cat$beta_pvalue), by = "name") %>%
  transmute(name, 
            value = round(exp(value.x), 2), 
            pvalue = round(value.y, 2))
}
