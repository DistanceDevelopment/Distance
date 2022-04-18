# process the ddf object
# other TODO:
#   - handle grouped/ungrouped estimation (all one or the other??)
#   - fix Total degrees of freedom
#   - left truncation
#   - dummy ddf
dht2_process_ddf <- function(ddf, convert_units, er_est, strat_formula){

  # if we don't have a list, make a list
  if(any(class(ddf) != "list")){
    ddf <- list(ddf)
  }

  # we can have a different unit conversion per detection function, as a treat
  if(length(convert_units) != 1){
    if(length(convert_units)!=length(ddf)){
      stop("convert_units must be either 1 number or have as many entries as there are detection functions")
    }else{
      convert_units <- rep(convert_units, length(ddf))
    }
  }

  # we can have a different ER estimators per detection function, as a treat
  # only check this if the par was set, else defaults get used
  if(!attr(er_est, "missing")){
    if(length(er_est) == 1){
      er_est <- rep(er_est, length(ddf))
      attr(er_est, "missing") <- FALSE
    }else if(length(er_est)!=length(ddf)){
      stop("er_est must be either 1 number or have as many entries as there are detection functions")
    }
  }

  # storage for the "distance" data
  bigdat <- c()
  obj_keep <- c()
  # storage for summaries
  ddf_summary <- list()
  transect_data <- data.frame(ddf_id = 1:length(ddf),
                              transect_type = rep(NA, length(ddf)),
                              er_est = rep(NA, length(ddf)),
                              df_width = rep(NA, length(ddf)))

  # just bad vibes below this point...
  for(i in seq_along(ddf)){

    this_ddf <- ddf[[i]]
    # just get the ds model if we have Distance::ds output
    if(inherits(this_ddf, "dsmodel")){
      this_ddf <- this_ddf$ddf
    }

    # drop unused levels of factors
    this_ddf$data <- droplevels(this_ddf$data)

    # only keep observations within the truncation
    obj_keep <- c(obj_keep, this_ddf$data$object[this_ddf$data$distance <=
                                                  this_ddf$meta.data$width &
                                                 this_ddf$data$distance >=
                                                  this_ddf$meta.data$left])
    this_bigdat <- this_ddf$data[this_ddf$data$object %in% obj_keep, ]

    # get probabilities of detection
    this_bigdat$p <- predict(this_ddf)$fitted

    # ensure as.factor in formula are propagated to the data
    this_bigdat <- safe_factorize(strat_formula, this_bigdat)

    # get variance estimation
    if(attr(er_est, "missing")){
      er_estl <- er_est[as.numeric(this_ddf$meta.data$point)+1]
    }else{
      er_estl <- er_est[i]
    }
    # transect data
    transect_data[i,] <- data.frame(ddf_id = i,
                                    # transect type
                                    transect_type=if(this_ddf$meta.data$point){
                                                      "point"}else{"line"},
                                    # ER variance estimation
                                    er_est = er_estl,
                                    # apply unit conversion to truncations
                                    df_width = this_ddf$meta.data$width *
                                               convert_units)

    # add a detection function identifier for this bit of the data
    this_bigdat$ddf_id <- i

    # put that back
    ddf[[i]] <- this_ddf
    bigdat <- rbind.data.frame(bigdat, this_bigdat)
  }

  # total number of data used to fit the detection functions
  transect_data$n_ddf <- sum(unlist(lapply(ddf, function(x) length(x$fitted))))
  # total number of parameters in the detection function
  transect_data$n_par <- sum(unlist(lapply(ddf, function(x) length(x$par))))

  if(any(table(bigdat$object) > 1)){
    stop("object column but be unique over all data")
  }

  list(ddf = ddf,
       bigdat = bigdat,
       obj_keep = obj_keep,
       transect_data = transect_data,
       summary  = ddf_summary)
}
