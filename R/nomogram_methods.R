#' Nomogram methods for cureit objects
#'
#' @param survival Logical indicating whether or not to output the nomogram
#' based on the survival submodel. Defaults to `TRUE`.
#' @param cure Logical indicating whether or not to output the nomogram
#' based on the cured probability submodel. Defaults to `TRUE`.

#' @param time Numeric vector of times to obtain survival probability estimates at

#'
#' @name nomogram_methods_cureit
#' @return a tibble
#' @family cureit() functions
#' @examples
#' c <- cureit(surv_formula = Surv(ttdeath, death) ~ age + grade, 
#' cure_formula = ~ age + grade,  data = trial)
#'

#' nomogram.cureit(x = c,time=300)

NULL

# nomogram
#' @rdname nomogram_methods_cureit
#' 
#' @export
#' @family cureit
nomogram.cureit <- function(x,
                            survival = TRUE, 
                            cure = TRUE,
                            time = NULL,
                            angle = 0,
                            ...) {
  
  # Data checks -------
  if (is.null(time)) {
    stop("Must specify at least one time point via `time=`.", call. = FALSE)
  }
  if (!is.null(time) && any(time < 0)) {
    stop("`time=` must be non-negative.", call. = FALSE)
  }
  if (length(time) > 1){
    stop("`time=` cannot be a vector.", call. = FALSE)
  }
  
  # get processed data & linear predictor
  processed <- cureit_mold(x$surv_formula, x$cure_formula, x$data)
  predicted_lp <- predict(x, times=time,method="lp")
  
  if (survival){
    
    surv.covnames <- all.vars(x$surv_formula[[3]])
    surv.factors <- names(x$surv_xlevels)
    surv.continuous <- setdiff(surv.covnames,surv.factors)
    num.levels <- sapply(x$surv_xlevels,length)
    
    # both continuous and categorical 
    if (length(surv.factors) > 0 & length(surv.continuous) > 0){
      surv_nmmat <- data.frame(variables=c(unlist(mapply(rep,surv.factors,num.levels)),
                                           unlist(mapply(rep,surv.continuous,10))),
                               levels=c(unlist(x$surv_xlevels),rep(NA,length(surv.continuous)*10)),
                               lp=0,
                               lpscale=0,
                               type=c(unlist(mapply(rep,rep("factor",length(surv.factors)),num.levels)),
                                      unlist(mapply(rep,rep("continuous",length(surv.continuous)),10)))
      )
      
    # categorical only 
    }else if(length(surv.factors) > 0 & length(surv.continuous) == 0){
      surv_nmmat <- data.frame(variables=c(unlist(mapply(rep,surv.factors,num.levels))),
                               levels=c(unlist(x$surv_xlevels)),
                               lp=0,
                               lpscale=0,
                               type=c(unlist(mapply(rep,rep("factor",length(surv.factors)),num.levels)))
      )
    # continuous only 
    }else if (length(surv.factors) == 0 & length(surv.continuous) > 0){
      surv_nmmat <- data.frame(variables=c(unlist(mapply(rep,surv.continuous,10))),
                               levels=c(rep(NA,length(surv.continuous)*10)),
                               lp=0,
                               lpscale=0,
                               type=c(unlist(mapply(rep,rep("continuous",length(surv.continuous)),10)))
      )
    }
    
    surv_nmmat$levels <- gsub(" ", "_", surv_nmmat$levels)
    surv_nmmat$combined <- ifelse(is.na(surv_nmmat$levels),
                                  surv_nmmat$variables,
                                  janitor::make_clean_names(
                                    paste0(surv_nmmat$variables,surv_nmmat$levels)))
    
    # scaling constant -------
    covmin <- apply(processed$surv_processed$predictors,2,min, na.rm = TRUE)
    covmax <- apply(processed$surv_processed$predictors,2,max, na.rm = TRUE)
    
    # variable range * survival coefficients 
    lpdiff <- (covmax-covmin)*x$surv_coefs
    
    # largest lp value needed among all
    nm_scale <- 100/max(abs(lpdiff))
    lpscale <- lpdiff * nm_scale
    
    
    for (i in 1:length(lpdiff)){
      type <- unique(surv_nmmat$type[surv_nmmat$combined == names(lpdiff)[i]])
      if (type=="factor"){
        surv_nmmat$lp[surv_nmmat$combined == names(lpdiff)[i]] = lpdiff[i]
        surv_nmmat$lpscale[surv_nmmat$combined == names(lpdiff)[i]] = lpscale[i]
        
      # continuous- create breakpoints 
      }else if(type=="continuous"){
        surv_nmmat$lp[surv_nmmat$combined == names(lpdiff)[i]] = lpdiff[i]
        surv_nmmat$lpscale[surv_nmmat$combined == names(lpdiff)[i]] = seq(0,lpscale[i],length.out=10)
        surv_nmmat$levels[surv_nmmat$combined == names(lpdiff)[i]] = pretty(covmin[names(lpdiff)[i]]:covmax[names(lpdiff)[i]],10)
      }
    }
    surv_nmmat$levels <- gsub("_", " ", surv_nmmat$levels)
    
    for (var in unique(surv_nmmat$variables)){
      if (min(surv_nmmat$lpscale[surv_nmmat$variables==var]) < 0){
        surv_nmmat$lpscale[surv_nmmat$variables==var] = surv_nmmat$lpscale[surv_nmmat$variables==var] - min(surv_nmmat$lpscale[surv_nmmat$variables==var])
      }
    }
    
    surv.lp <- pretty(range(predicted_lp$lp_surv_model, na.rm = TRUE), n = 10)
    surv.lp.scaled <- (surv.lp-min(surv.lp)) * nm_scale
    
    points = pretty(0:100, 10)
    df_points <- as.data.frame(points) %>%
      transmute(x = points,
                levels = as.character(points), 
                y = "Points")
    
    upper_range <-  max(surv.lp.scaled, 100, na.rm = TRUE)
    total_points <- pretty(c(0, upper_range), n=10)
    upper_range_pretty <- max(total_points, na.rm = TRUE)
    
    df_lp_surv <- as.data.frame(surv.lp) %>%
      transmute(x = surv.lp.scaled, 
                levels = as.character(surv.lp), 
                y = "Linear Predictor") %>%
      mutate(x = (x/upper_range_pretty)*100)
    
    df_tp_surv <- as.data.frame(total_points) %>%
      mutate(levels = as.character(round(total_points, 0)), 
             x = scales::rescale(total_points, to = c(0, 100)), 
             y = "Uncured survival: \nTotal Points") %>%
      select(-total_points)
    
    baseline <- cbind.data.frame(s = x$smcure$s, times = x$smcure$Time) %>%
      mutate(diff = abs(times-time)) %>%
      filter(diff == min(abs(diff)))
    s0 <- baseline$s[1]
    
    df_pred_surv <- as.data.frame(surv.lp.scaled) %>%
      transmute(x = (surv.lp.scaled/upper_range_pretty)*100, 
                levels = as.character(round(s0^exp(surv.lp), 1)), 
                y = paste0("Uncured Survival: \nProbability, T=",time))
    
    df_cov_surv <- data.frame(x=surv_nmmat$lpscale,levels=surv_nmmat$levels,y=paste0("Uncured Survival: \n",
                                                                                surv_nmmat$variables))
    
    all_surv <- bind_rows(df_cov_surv, df_tp_surv, df_pred_surv)
    all_surv$model <- "surv"

  }
  
  if (cure){
    
    cure.covnames <- all.vars(x$cure_formula[[2]])

    cure.factors <- names(x$cure_xlevels)
    cure.continuous <- setdiff(cure.covnames,cure.factors)
    num.levels <- sapply(x$cure_xlevels,length)
    if (length(cure.factors) > 0 & length(cure.continuous) > 0){
      cure_nmmat <- data.frame(variables=c(unlist(mapply(rep,cure.factors,num.levels)),
                                           unlist(mapply(rep,cure.continuous,10))),
                               levels=c(unlist(x$cure_xlevels),rep(NA,length(cure.continuous)*10)),
                               lp=0,
                               lpscale=0,
                               type=c(unlist(mapply(rep,rep("factor",length(cure.factors)),num.levels)),
                                      unlist(mapply(rep,rep("continuous",length(cure.continuous)),10)))
      )
    }else if(length(cure.factors) > 0 & length(cure.continuous) == 0){
      cure_nmmat <- data.frame(variables=c(unlist(mapply(rep,cure.factors,num.levels))),
                               levels=c(unlist(x$cure_xlevels)),
                               lp=0,
                               lpscale=0,
                               type=c(unlist(mapply(rep,rep("factor",length(cure.factors)),num.levels)))
      )
    }else if (length(cure.factors) == 0 & length(cure.continuous) > 0){
      cure_nmmat <- data.frame(variables=c(unlist(mapply(rep,cure.continuous,10))),
                               levels=c(rep(NA,length(cure.continuous)*10)),
                               lp=0,
                               lpscale=0,
                               type=c(unlist(mapply(rep,rep("continuous",length(cure.continuous)),10)))
      )
    }
    cure_nmmat$levels <- gsub(" ", "_", cure_nmmat$levels)
    cure_nmmat$combined <- ifelse(is.na(cure_nmmat$levels),
                                  cure_nmmat$variables,
                                  janitor::make_clean_names(paste0(cure_nmmat$variables,cure_nmmat$levels)))
    
    covmin <- apply(processed$cure_processed$predictors,2,min, na.rm = TRUE)
    covmax <- apply(processed$cure_processed$predictors,2,max, na.rm = TRUE)
    lpdiff <- (covmax-covmin)*x$cure_coefs[-1]
    nm_scale <- 100/max(abs(lpdiff))
    lpscale <- lpdiff * nm_scale
    
    for (i in 1:length(lpdiff)){
      type <- unique(cure_nmmat$type[cure_nmmat$combined == names(lpdiff)[i]])
      if (type=="factor"){
        cure_nmmat$lp[cure_nmmat$combined == names(lpdiff)[i]] = lpdiff[i]
        cure_nmmat$lpscale[cure_nmmat$combined == names(lpdiff)[i]] = lpscale[i]
      }else if(type=="continuous"){
        cure_nmmat$lp[cure_nmmat$combined == names(lpdiff)[i]] = lpdiff[i]
        cure_nmmat$lpscale[cure_nmmat$combined == names(lpdiff)[i]] = seq(0,lpscale[i],length.out=10)
        cure_nmmat$levels[cure_nmmat$combined == names(lpdiff)[i]] = pretty(covmin[names(lpdiff)[i]]:covmax[names(lpdiff)[i]],10)
      }
    }
    cure_nmmat$levels <- gsub("_", " ", cure_nmmat$levels)
    
    for (var in unique(cure_nmmat$variables)){
      if (min(cure_nmmat$lpscale[cure_nmmat$variables==var]) < 0){
        cure_nmmat$lpscale[cure_nmmat$variables==var] = cure_nmmat$lpscale[cure_nmmat$variables==var] - min(cure_nmmat$lpscale[cure_nmmat$variables==var])
      }
    }
    
    cure.lp <- pretty(range(predicted_lp$lp_cure_model, na.rm = TRUE), n = 10)
    cure.lp.scaled <- (cure.lp-min(cure.lp, na.rm = TRUE)) * nm_scale
    
    # points = pretty(0:100, 10)
    # df_points <- as.data.frame(points) %>%
    #   transmute(x = points,
    #             levels = as.character(points), 
    #             y = "Cured probability: \nPoints")
    
    upper_range <-  max(cure.lp.scaled, 100, na.rm = TRUE)
    total_points <- pretty(c(0, upper_range), n=10)
    upper_range_pretty <- max(total_points, na.rm = TRUE)
    
    df_lp_cure <- as.data.frame(cure.lp) %>%
      transmute(x = cure.lp.scaled, 
                levels = as.character(cure.lp), 
                y = "Linear Predictor") %>%
      mutate(x = (x/upper_range_pretty)*100)
    
    df_tp_cure <- as.data.frame(total_points) %>%
      mutate(levels = as.character(round(total_points, 0)), 
             x = scales::rescale(total_points, to = c(0, 100)), 
             y = "Cured probability: \nTotal Points") %>%
      select(-total_points)
    
    df_pred_cure <- as.data.frame(cure.lp.scaled) %>%
      transmute(x = (cure.lp.scaled/upper_range_pretty)*100, 
                levels = as.character(round(1-exp(cure.lp)/(1+exp(cure.lp)), 1)), 
                y = "Cured probability")
    
    df_cov_cure <- data.frame(x=cure_nmmat$lpscale,levels=cure_nmmat$levels,y=paste0("Cured probability: \n",
                                                                                cure_nmmat$variables))
    
    all_cure <- bind_rows(df_cov_cure, df_tp_cure, df_pred_cure)
    all_cure$model <- "cure"
    
  }
  
  df_points$model <- "NA"
  all <- bind_rows(df_points, all_cure, all_surv)
  
  all <- all %>%
    mutate(y = forcats::fct_relevel(y, unique(all$y))) %>%
    mutate(y = forcats::fct_rev(y)) %>% 
    mutate(model = forcats::fct_relevel(model, unique(all$model)))
  
  p1 <- all %>%
    ggplot(aes(x = x, y = y)) + geom_line(aes(color=model)) +
    geom_point(aes(color=model)) + 
    geom_text(aes(label = levels), vjust = 1.5, angle=angle)  + ylab(" ") + xlab(" ") +
    # ggtitle("Estimated cureival for Uncured") +
    theme_minimal() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank()) +
    scale_color_manual(values=c("black","red","blue"))+ 
    guides(color="none")
  
  p1

}