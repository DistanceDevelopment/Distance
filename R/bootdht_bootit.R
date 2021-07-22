# function to do a single bootstrap iteration
bootit <- function(bootdat, models, our_resamples, summary_fun,
                   convert.units, pb, multipliers_fun, sample_label,
                   select_adjustments, ...){

  # get resampled data
  bootdat <- bootdht_resample_data(bootdat, our_resamples)

  aics <- rep(NA, length(models))
  for(i in seq_along(models)){
    model <- models[[i]]

    # setup the call to ds
    df_call <- model$call

    # don't calculate standard errors
    df_call$dht.se <- FALSE

    # if we want the number of adjustments to be selected each iteration...
    if(select_adjustments){
      df_call$order <- NULL
    }
    # insert the new data into the model
    df_call$data <- bootdat
    if(!is.null(convert.units)){
      df_call$convert.units <- convert.units
    }

    # fit that and update what's in models
    models[[i]] <- try(suppressMessages(eval(df_call)),
                       silent=TRUE)

    if(any(class(models[[i]]) == "try-error") ||
       is.null(models[[i]]$dht)){
      # if the model failed, return NA
      aics[i] <- NA
    }else{
      # if that wasn't bad, grab the AIC
      aics[i] <- AIC(models[[i]])$AIC

      # apply multipliers
      rates <- unique(bootdat[, c("Region.Label", "rate")])
      # if there's only one stratum, it's always named Total
      if(nrow(rates)==1) rates$Region.Label <- "Total"
      # take the product of all the functional multipliers and the
      # non-functional ones too
      rates$rate <- rates$rate *
                      prod(unlist(lapply(multipliers_fun, function(f) f())))

      # now need to merge the rates object onto the results, and re-scale
      # abundance/density estimates for individuals
      indN <- merge(models[[i]]$dht$individuals$N, rates,
                    by.x="Label", by.y="Region.Label",
                    all.x=TRUE)
      indD <- merge(models[[i]]$dht$individuals$D, rates,
                    by.x="Label", by.y="Region.Label",
                    all.x=TRUE)
      if(nrow(indN)==1){
        models[[i]]$dht$individuals$N$Estimate <- indN$Estimate/indN$rate
        models[[i]]$dht$individuals$D$Estimate <- indD$Estimate/indD$rate
      }else{
        nN <- length(indN$Estimate)-1
        models[[i]]$dht$individuals$N$Estimate <- c(indN$Estimate[1:nN]/
                                                    indN$rate[1:nN],
                                                    sum(indN$Estimate[1:nN]/
                                                        indN$rate[1:nN]))
        models[[i]]$dht$individuals$D$Estimate <- c(indD$Estimate[1:nN]/
                                                    indD$rate[1:nN],
                                                    sum(indD$Estimate[1:nN]/
                                                        indD$rate[1:nN]))
      }
      # ... and for clusters
      if(any(names(models[[i]]$dht)=="clusters")){
        clN <- merge(models[[i]]$dht$clusters$N, rates, by="Region.Label",
                     all.x=TRUE)
        clD <- merge(models[[i]]$dht$clusters$D, rates, by="Region.Label",
                     all.x=TRUE)
        if(nrow(indN)==1){
          models[[i]]$dht$clusters$N$Estimate <- clN$Estimate/clN$rate
          models[[i]]$dht$clusters$D$Estimate <- clD$Estimate/clD$rate
        }else{
          models[[i]]$dht$clusters$N$Estimate <- c(clN$Estimate[1:nN]/
                                                      clN$rate[1:nN],
                                                      sum(clN$Estimate[1:nN]/
                                                          clN$rate[1:nN]))
          models[[i]]$dht$clusters$D$Estimate <- c(clD$Estimate[1:nN]/
                                                      clD$rate[1:nN],
                                                      sum(clD$Estimate[1:nN]/
                                                          clD$rate[1:nN]))
        }
      }
    }
  }

  # update progress bar
  pb$increment(pb$pb)

  if(all(is.na(aics))){
    # if no models fitted, return NA
    ret <- NA
    class(ret) <- "bootstrap_failure"
    return(ret)
  }else{
    fit <- models[[which.min(aics)]]
    # handle errors
    if(any(class(fit) == "try-error") ||
       any(is.na(fit$ddf$hessian))){
      ret <- NA
      class(ret) <- "bootstrap_failure"
      return(ret)
    }else{
      return(summary_fun(fit$dht, fit$ddf))
    }
  }
}

