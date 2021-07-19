# function to do a single bootstrap iteration
bootit <- function(bootdat, models, our_resamples, summary_fun,
                   convert.units, pb, multipliers_fun, sample_label,
                   select_adjustments, ...){
  # sample at the right levels
  for(sample_thingo in our_resamples){
    # what are the possible samples at this level
    levs <- unique(bootdat[[sample_thingo]])
    nlevs <- length(levs)
    levs <- sample(levs, nlevs, replace=TRUE)

    # make a new data frame with the correct number of replicates of the
    # per-stratum data in it
    bootdat <- lapply(levs, function(x){
      bootdat[bootdat[[sample_thingo]] == x, ]
    })
    # make a special index to make unique IDs later
    iind <- rep(seq_len(length(bootdat)), lapply(bootdat, nrow))
    # make list of data.frames into one frame
    bootdat <- do.call("rbind", bootdat)
    # put that ID in there
    bootdat[[paste0(sample_thingo, "_ID")]] <- iind
  }

  # need unique object IDs
  bootdat$object <- seq_len(nrow(bootdat))
  # get the sample labels right
  bootdat$Sample.Label <- paste0(bootdat[[sample_label]], "-",
                                 bootdat[[paste0(sample_thingo, "_ID")]])

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

    if(any(class(models[[i]]) == "try-error")){
      # if the model failed, return NA
      aics[i] <- NA
    }else{
      # if that wasn't bad, grab the AIC
      aics[i] <- AIC(models[[i]])$AIC

      # apply multipliers
      rates <- unique(bootdat[, c("Region.Label", "rate")])
      # if there's only one stratum, it's always named Total
      if(nrow(rates)==1) rates$Region.Label <- "Total"
      # take the product of all the functional multipliers
      rates$rate <-  prod(unlist(lapply(multipliers_fun, function(f) f())))

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
        models[[i]]$dht$individuals$N$Estimate <- c(indN$Estimate/indN$rate,
                                           sum(indN$Estimate/indN$rate))
        models[[i]]$dht$individuals$D$Estimate <- c(indD$Estimate/indD$rate,
                                           sum(indD$Estimate/indD$rate))
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
          models[[i]]$dht$clusters$N$Estimate <- c(clN$Estimate/clN$rate,
                                          sum(clN$Estimate/clN$rate))
          models[[i]]$dht$clusters$D$Estimate <- c(clD$Estimate/clD$rate,
                                          sum(clD$Estimate/clD$rate))
        }
      }
    }
  }

  # update progress bar
  pb$increment(pb$pb)

  if(all(is.na(aics))){
    # if no models fitted, return NA
    nbootfail <<- nbootfail + 1
    return(NA)
  }else{
    fit <- models[[which.min(aics)]]
    # handle errors
    if(any(class(fit) == "try-error") ||
       any(is.na(fit$ddf$hessian))){
      nbootfail <<- nbootfail + 1
      return(NA)
    }else{
      return(summary_fun(fit$dht, fit$ddf))
    }
  }
}
