#' Handle multipliers for abundance estimation
#'
#' This is an internal function.
#' @noRd
dht2_multipliers <- function(multipliers, bigdat){

  # merge-in multipliers
  if(!is.null(multipliers)){
    if(!is.list(multipliers)){
      stop("multipliers must be a list")
    }
    if(!all(names(multipliers) %in% c("creation", "decay")) |
       is.null(names(multipliers))){
      stop("Multipliers must be named \"creation\" and \"decay\"")
    }
    if(length(multipliers)>2){
      stop("Only one creation and one decay rate may be provided")
    }

    # if only creation or decay were included, add a dummy version to
    # make things easier below
    if(length(multipliers)==1){
      miss_mult <- setdiff(c("decay", "creation"), names(multipliers))
      multipliers[[miss_mult]] <- data.frame(rate = 1,
                                             df   = 0,
                                             SE   = 0,
                                             CV   = 0)
    }

    for(ii in names(multipliers)){
      # check that length is correct
      if(nrow(multipliers[[ii]])>1 &
         (all(names(multipliers[[ii]]) %in% c("rate", "df", "SE")) |
          !any(names(multipliers[[ii]]) %in% names(bigdat)))){
        stop("Multirow multipliers need column to link to the data")
      }
      if(!is.data.frame(multipliers[[ii]])){
        stop(paste0("multipliers[[", ii, "]] must be a data.frame"))
      }
      # check multipliers has at least a rate column
      if(!("rate" %in% names(multipliers[[ii]]))){
        stop(paste("You need at least a column named \"rate\" in",
                   names(multipliers)[ii], "multiplier"))
      }

      if(is.null(multipliers[[ii]]$df)){
        # this deals with the no df case
        multipliers[[ii]]$df <- Inf
      }
      if(is.null(multipliers[[ii]]$SE)){
        multipliers[[ii]]$SE <- 0
      }

      # calculate CV
      multipliers[[ii]]$CV <- (multipliers[[ii]]$SE/multipliers[[ii]]$rate)

      if(ii=="decay" && multipliers[[ii]]$SE!=0){
        # since we are dividing, use the sandwich estimator,
        # var(1/x) = 1/x^2 * var(x) * 1/x^2 => se(1/x) = se(x)/x^2
        multipliers[[ii]]$SE <- multipliers[[ii]]$SE/multipliers[[ii]]$rate^2
      }
      for(iii in c("rate", "SE", "CV", "df")){
        multipliers[[ii]][[paste0(ii, "_", iii)]] <- multipliers[[ii]][[iii]]
        multipliers[[ii]][[iii]] <- NULL
      }

      bigdat <- merge(bigdat, multipliers[[ii]], all.x=TRUE)
    }
    # get the "final" estimates
    bigdat$rate <- bigdat$creation_rate*bigdat$decay_rate
    bigdat$rate_SE <- sqrt(bigdat$creation_SE^2 + bigdat$decay_SE^2)
    bigdat$rate_df <- bigdat$creation_df +  bigdat$decay_df
    bigdat$rate_CV <- sqrt(bigdat$creation_CV^2 + bigdat$decay_CV^2)

    # remove the creation/decay columns that we created
    del_names <- c(names(multipliers[["decay"]]),
                   names(multipliers[["creation"]]))
    del_names <- del_names[(grepl("^creation_", del_names) |
                             grepl("^decay_", del_names))]
    bigdat[del_names] <- NULL
    attr(bigdat, "multipliers") <- TRUE
  }else{
    # setup "fake" data for when we don't have multipliers
    # this makes the calculations cleaner below
    bigdat <- bigdat %>%
      mutate(Nc_cuecorrected = NA,
             rate = 1,
             rate_df = 1,
             rate_SE = 0,
             rate_CV = 0)
    attr(bigdat, "multipliers") <- FALSE
  }

  # we use var everywhere else!
  bigdat$rate_var <- bigdat$rate_SE^2
  bigdat$rate_SE <- NULL

  return(bigdat)
}
