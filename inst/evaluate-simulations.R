source(here::here("R/cv.cureitlasso.R"))
source(here::here("R/cureitlasso.R"))
source(here::here("R/simulasso.R"))
source(here::here("R/coxsplit.R"))

# Calculate Performance Metrics on Simulations ----------------------------------

# * Load all relevant data -------
load(here::here("inst", "cv_fits_50.RData"))
load(here::here("inst", "cv_fits_100.RData"))
load(here::here("inst", "sim_data_list.RData"))
load(here::here("inst", "sim_valid_data_list.RData"))

# * Calculate Sensitivity of Variable Selection in Training Fits -------
x <- for (i in 1:50) {
  fit <- cv_fits_list[[i]]
  #fit <- x$fit
  dat <- sim_data_list[[i]]
  
  sen_cox_min <- mean(dat$id_cox %in% which(fit$fit[[fit$index$min[1]]][[fit$index$min[2]]]$beta!=0))
  print(sen_cox_min)
  
}


hist(x)

# Sensitivity and specificity (for both min and 1se choices of parameters) ---

sen_cox_min <- mean(dat$id_cox %in% which(fit$fit[[fit$index$min[1]]][[fit$index$min[2]]]$beta!=0))


spe_cure_min <- 1-mean(setdiff(1:200,dat$id_cure) %in% which(fit$fit[[fit$index$min[1]]][[fit$index$min[2]]]$alpha!=0))
spe_cox_min <- 1-mean(setdiff(1:200,dat$id_cure) %in% which(fit$fit[[fit$index$min[1]]][[fit$index$min[2]]]$beta!=0))

sen_cure_1se <- mean(dat$id_cure %in% which(fit$fit[[fit$index$`1se`[1]]][[fit$index$`1se`[2]]]$alpha!=0))
sen_cox_1se <- mean(dat$id_cox %in% which(fit$fit[[fit$index$`1se`[1]]][[fit$index$`1se`[2]]]$beta!=0))
spe_cure_1se <- 1-mean(setdiff(1:200,dat$id_cure) %in% which(fit$fit[[fit$index$`1se`[1]]][[fit$index$`1se`[2]]]$alpha!=0))
spe_cox_1se <- 1-mean(setdiff(1:200,dat$id_cure) %in% which(fit$fit[[fit$index$`1se`[1]]][[fit$index$`1se`[2]]]$beta!=0))




# * Get it to work with yardstick -----
library(yardstick)

vec <- c(TRUE, TRUE, FALSE)
df <- data.frame(
  pred = factor(as.integer(dat$id_cure %in% which(fit$fit[[fit$index$min[1]]][[fit$index$min[2]]]$alpha!=0)))
  # Presence = as.integer(vec)  # Convert logical values to integers (TRUE -> 1, FALSE -> 0)
) %>%
  mutate(obs = factor(rep(1, length = nrow(.)), levels = c("0", "1")))

sens(df, obs, pred, event_level = "second")
sen_cure_min <- mean(dat$id_cure %in% which(fit$fit[[fit$index$min[1]]][[fit$index$min[2]]]$alpha!=0))



# * Brier Score -------------------------------------------------------------
#' - `ti`- simulated survival times
#' - `d`- simulated event indicator
#' - `x` simulated data n x p matrix or data.frame
#' - `id_cure` - IDs of selected cure parameters
#' - `id_cure` - IDs of selected cox parameters
#' - `coefs_cure` - Coefficients of cure parameters
#' - `coefs_cox` - Coefficients of cox parameters
#' - `tj` - is unique event times sorted


for (i in 1:length(sim_data_list)) {
  
  dat <- sim_data_list[[i]]
  dat_valid = sim_valid_data_list[[i]]
  fit <- cv_fits_list[[i]]
  
  # for training set
  ti <- dat$t
  di <- dat$d
  xi <- dat$x
  # unique event times sorted
  tj <- fit$fit[[fit$index$`min`[1]]][[fit$index$`min`[2]]]$tj
  
  # for test set
  ti_valid <- dat_valid$t
  di_valid <- dat_valid$d
  xi_valid <- dat_valid$x
  tj_valid <- sort(ti_valid[di_valid==1])
  
  
  # get final fit (hyperparam with min error) for cure model (fits used train data)
  final_select_model = fit$fit[[fit$index$`min`[1]]][[fit$index$`min`[2]]]
  fitcure <- final_select_model$fitcure
  
  # use this model to make predictions on data (fitted relative-risk for "cox)
  predcure <- predict(fitcure, newx=xi, s=min(fitcure$lambda), type="response")
  predcure_valid <- predict(fitcure, newx = xi_valid, s=min(fitcure$lambda), type="response")
  
  # get final fit of min coxcure model
  fitcox <- final_select_model$fitcox
  
  # Get Data from Selected "Final" Prediction Model -----
  # use this model to make predictions on data
  
  # get linear predictors from cox
  predsurvexp <- predict(fitcox,
                         newx=xi,
                         s = min(fitcox$lambda),
                         type="response")
  
  predsurvexp_valid <- predict(fitcox,
                               newx=xi_valid,
                               s = min(fitcox$lambda),
                               type="response")
  
  haz <- final_select_model$haz
  cumhaz <- final_select_model$cumhaz
  
  # Get Kaplan Meier estimate of actual data
  predcens <- survfit(Surv(ti,1-di)~1)$surv
  tcens <- survfit(Surv(ti,1-di)~1)$time
  # 
  # predcens <- survfit(Surv(ti,1-di)~1)$surv
  # tcens <- survfit(Surv(ti,1-di)~1)$time
  
  # times at which to get score
  tbrier <- sort(ti)
  brier <- rep(NA,length(tbrier))
  
  tbrier_valid <- sort(ti_valid)
  brier_valid <- rep(NA,length(tbrier_valid))
  
  # Train Data----
  for (l in 1:length(tbrier)){
    
    predsurv <- rep(NA,length(ti)) # Conditional survival prob for all patients at tbrier[l]
    
    if (tbrier[l] >= min(tj)){
      
      ids <- which(tj == max(tj[tj <= tbrier[l]]))
      
      # get predicted survival probability
      predsurv <- exp(-cumhaz[ids]*predsurvexp)
      
      # weights
      ipw <- as.numeric(ti > tbrier[l])*predcens[l] + as.numeric(ti <= tbrier[l])*predcens
      
    } else if (tbrier[l] < min(tj)){
      
      predsurv <- rep(1,length(ti))
      ipw <- 1
      
    }
    
    preds <- 1 - predcure + predcure*predsurv
    brier[l] <- mean(I(di==0)*(1 - preds)^2/(ipw+0.001) + I(di==1)*(0 - preds)^2/(ipw+0.001))
    
  }
  
  # Test Data----
  for (l in 1:length(tbrier)){
    
    predsurv_valid <- rep(NA,length(ti_valid)) # Conditional survival prob for all patients at tbrier[l]
    
    if (tbrier_valid[l] >= min(tj_valid)){
      
      ids_valid <- which(tj_valid == max(tj_valid[tj_valid <= tbrier_valid[l]]))
      
      # get predicted survival probability
      predsurv_valid <- exp(-cumhaz[ids_valid]*predsurvexp_valid)
      
      # weights
      ipw_valid <- as.numeric(ti_valid > tbrier_valid[l])*predcens[l] + as.numeric(ti_valid <= tbrier_valid[l])*predcens
      
    }else if (tbrier_valid[l] < min(tj_valid)){
      
      predsurv_valid <- rep(1,length(ti_valid))
      ipw_valid <- 1
      
    }
    
    preds_valid <- 1 - predcure_valid + predcure_valid*predsurv_valid
    brier_valid[l] <- mean(I(di_valid==0)*(1 - preds_valie)^2/(ipw_valid + 0.001) + I(di_valid==1)*(0 - preds_valid^2/(ipw_valid+0.001))
                           
                           
                           
                           all_brier = cbind.data.frame("tbrier" = tbrier, "brier" = brier, "brier_valid" = brier_valid)
                           return(all_brier)
                           # plot(tbrier,brier,type="S")
                           
                           
  }
}
