library(cureit)
library(ggsurvfit)

load(here::here("inst", "deid_data_2022-09-02[1526].RData"))

# Summary ------
fit <- survfit(Surv(os_months, dod) ~ 1, data = df_sel)

ggsurvfit(fit) + ylim(0, 1) + 
  add_risktable() + theme_classic()


# Test Inference
df <- df_sel %>% select(os_months,
                            dod, variant_name_specify) %>%
  na.omit()

p <- cureit(surv_formula = Surv(os_months, dod) ~ variant_name_specify,
            cure_formula = ~ variant_name_specify,
            data = df,
            nboot = 2)

times <- seq(5, 24, 0.5)
bootbrier <- Brier_inference_bootstrap(p, times = times, nboot = 2)
plot(times, bootbrier$brier, type = "s", xlab = "Time", ylab = "Brier score")
lines(times, bootbrier$brier_2.5, type = "s", lty = 2)
lines(times, bootbrier$brier_97.5, type = "s", lty = 2)
lines(times, bootbrier$brier_cox, type = "s", col = "red")
lines(times, bootbrier$brier_cox_2.5, type = "s", lty = 2, col = "red")
lines(times, bootbrier$brier_cox_97.5, type = "s", lty = 2, col = "red")
legend("topleft", c("Cure model", "Cox model"), col = c("black", "red"), lty = 1)


# variant name specify -----
set.seed(103022)
set.seed(10)
x <- cureit(surv_formula = Surv(os_months, dod) ~ variant_name_specify,
            cure_formula = ~ variant_name_specify, data = df_sel,
            nboot = 200,
            eps = 0.05)

nomogram(x,time=400)
pred <- predict(x,method="prob",times=seq(10,390,10),brier=TRUE,cox=TRUE)

plot(seq(10,390,10),pred$brier,type="s",xlab="Time",ylab="Brier score",ylim=c(0,0.15))
lines(seq(10,390,10),pred$brier_cox,type="s",col="red")
legend("topleft",c("Cure model","Cox model"),col=c("black","red"),lty=1)

pred1 <- pred

# Multivariate:  Variant/Size------

set.seed(103022)

df <- df_sel %>%
  select(os_months, dod, variant_name_specify, max_tumor_size) %>%
  na.omit()

x <- cureit(surv_formula = Surv(os_months, dod) ~ variant_name_specify + log(max_tumor_size + 1),
            cure_formula = ~ variant_name_specify+ log(max_tumor_size + 1),
            data = df,
            nboot = 100,
            eps = 0.05)

#nomogram(x,time=400)
pred <- predict(x,method="prob",
                times=seq(10,390,10),
                brier=TRUE,cox=TRUE)

plot(seq(10,390,10),pred$brier,type="s",xlab="Time",ylab="Brier score",ylim=c(0,0.15))
lines(seq(10,390,10),pred$brier_cox,type="s",col="red")
legend("topleft",c("Cure model","Cox model"),col=c("black","red"),lty=1)

pred2 <- pred

# Multivariate:  Variant/Size/Depth------

set.seed(103022)

df <- df_sel %>%
  select(os_months, dod, variant_name_specify, max_tumor_size, depth_dscrp) %>%
  na.omit()

x <- cureit(surv_formula = Surv(os_months, dod) ~ variant_name_specify + log(max_tumor_size + 1) + depth_dscrp,
            cure_formula = ~ variant_name_specify+ log(max_tumor_size + 1) + depth_dscrp,
            data = df,
            nboot = 100,
            eps = 0.05)

#nomogram(x,time=400)
pred <- predict(x,method="prob",
                times=seq(10,390,10),
                brier=TRUE,cox=TRUE)

plot(seq(10,390,10),pred$brier,type="s",xlab="Time",ylab="Brier score",ylim=c(0,0.15))
lines(seq(10,390,10),pred$brier_cox,type="s",col="red")
legend("topleft",c("Cure model","Cox model"),col=c("black","red"),lty=1)


pred3 <- pred

# Multivariate:  Variant/Size/Depth/Location ------

set.seed(103022)

df <- df_sel %>%
  select(os_months, dod, variant_name_specify, max_tumor_size, depth_dscrp, site1_name) %>%
  na.omit()

x <- cureit(surv_formula = Surv(os_months, dod) ~ variant_name_specify + log(max_tumor_size + 1) + depth_dscrp + site1_name,
            cure_formula = ~ variant_name_specify+ log(max_tumor_size + 1) + depth_dscrp + site1_name,
            data = df,
            nboot = 100,
            eps = 0.05)

#nomogram(x,time=400)
pred <- predict(x,method="prob",
                times=seq(10,390,10),
                brier=TRUE,cox=TRUE)

plot(seq(10,390,10),pred$brier,type="s",xlab="Time",ylab="Brier score",ylim=c(0,0.15))
lines(seq(10,390,10),pred$brier_cox,type="s",col="red")
legend("topleft",c("Cure model","Cox model"),col=c("black","red"),lty=1)

pred4 <- pred

# Compare ------

compare_brier <- function(pred) {
  pred$brier_cox - pred$brier
}


#all_pred <- list(pred1, pred2, pred3, pred4)
names(all_pred) <- c("variant", "variant_size", "variant_size_depth", "variant_size_depth_site")
brier_dif <- purrr::map_dfc(list(pred1, pred2, pred3, pred4), ~compare_brier(.x))

names(brier_dif) <- c("variant", "variant_size", "variant_size_depth",
                      "variant_size_depth_site")
brier_dif <- brier_dif %>%
  mutate(time =seq(10,390,10)) %>%
  tidyr::pivot_longer(-time)

library(ggplot2)
ggplot(brier_dif, aes(x = time, y = value, color = name)) + 
         geom_line() + theme_classic() + ylab("Brier Cox - Brier Cure")








