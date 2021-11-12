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
    df_call$dht_se <- FALSE

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
      which_ests <- c("N", "D")
      which_ests <- which_ests[which_ests %in%
                               names(models[[i]]$dht$individuals)]

      for(est_type in which_ests){
        ind <- merge(models[[i]]$dht$individuals[[est_type]], rates,
                     by.x="Label", by.y="Region.Label",
                     all.x=TRUE)
        if(nrow(ind)==1){
          models[[i]]$dht$individuals[[est_type]]$Estimate <- ind$Estimate/
                                                              ind$rate
        }else{
          nN <- length(ind$Estimate)-1
          models[[i]]$dht$individuals[[est_type]]$Estimate <-
            c(ind$Estimate[1:nN]/ind$rate[1:nN],
              sum(ind$Estimate[1:nN]/ind$rate[1:nN]))
        }
        # ... and for clusters
        if(any(names(models[[i]]$dht)=="clusters")){
          cl <- merge(models[[i]]$dht$clusters[[est_type]], rates,
                      by.x="Label", by.y="Region.Label",
                      all.x=TRUE)
          if(nrow(cl)==1){
            models[[i]]$dht$clusters[[est_type]]$Estimate <- cl$Estimate/
                                                             cl$rate
          }else{
            models[[i]]$dht$clusters[[est_type]]$Estimate <-
              c(cl$Estimate[1:nN]/cl$rate[1:nN],
                sum(cl$Estimate[1:nN]/cl$rate[1:nN]))
          }
        } # end clusters
      } # end loop over D, N (maybe)
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

