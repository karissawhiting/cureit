#### Evaluate the calibration by Brier score for the validation data
#### Line 79 - 123: Brier score for the cureit model

# ti <- dat_valid$t
# di <- dat_valid$d
# xi <- dat_valid$x
ti <- dat$t
di <- dat$d
xi <- dat$x
tj <- fit$fit[[fit$index$`min`[1]]][[fit$index$`min`[2]]]$tj

fitcure <- fit$fit[[fit$index$`min`[1]]][[fit$index$`min`[2]]]$fitcure
predcure <- predict(fitcure,newx=xi,s=min(fitcure$lambda),type="response")
fitcox <- fit$fit[[fit$index$`min`[1]]][[fit$index$`min`[2]]]$fitcox
predsurvexp <- predict(fitcox,newx=xi,s=min(fitcox$lambda),type="response")
haz <- fit$fit[[fit$index$`min`[1]]][[fit$index$`min`[2]]]$haz
cumhaz <- fit$fit[[fit$index$`min`[1]]][[fit$index$`min`[2]]]$cumhaz

predcens <- survfit(Surv(ti,1-di)~1)$surv
tcens <- survfit(Surv(ti,1-di)~1)$time

tbrier <- sort(ti)
brier <- rep(NA,length(tbrier))
ipw <- rep(NA,length(tbrier))
for (l in 1:length(tbrier)){

  predsurv <- rep(NA,length(ti)) # Conditional survival prob for all patients at tbrier[l]

  if (tbrier[l] >= min(tj)){

    ids <- which(tj == max(tj[tj <= tbrier[l]]))
    predsurv <- exp(-cumhaz[ids]*predsurvexp)
    ipw <- as.numeric(ti > tbrier[l])*predcens[l] + as.numeric(ti <= tbrier[l])*predcens

  }else if (tbrier[l] < min(tj)){

    predsurv <- rep(1,length(ti))
    ipw <- 1

  }
  
  ipw[l] <- ipw
  preds <- 1 - predcure + predcure*predsurv
  brier[l] <- mean(I(di==0)*(1 - preds)^2/(ipw+0.001) + I(di==1)*(0 - preds)^2/(ipw+0.001))

}

plot(tbrier,brier,type="S")