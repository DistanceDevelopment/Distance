bootdht_resample_data <- function(bootdat, our_resamples,
                                  stratum_label="Region.Label",
                                  sample_label="Sample.Label",
                                  obs_label="object"){

  # what are the unique combinations
  bootdat_samps <- unique(bootdat[, c(stratum_label, sample_label, obs_label)])

  # which strata?
  strata <- unique(bootdat_samps[[stratum_label]])
  # resample if that's what we're doing, else we keep all strata
  if(stratum_label %in% our_resamples){
    strata <- sample(strata, length(strata), replace=TRUE)
  }

  # which samples?
  # get all samples per stratum
  samps_per_strata <- by(bootdat_samps, bootdat_samps$Region.Label,
                         function(x) unique(x$Sample.Label))
  # restrict to those selected above (so there may be repeats of the
  # stratum name here, don't use names(samps_per_strata)!
  samps_per_strata <- samps_per_strata[strata]

  # resample if that's what we're doing, else we keep all samples
  samples <- unlist(samps_per_strata)
  if(sample_label %in% our_resamples){
    samples <- unlist(lapply(samps_per_strata, function(x){
      sample(x, length(x), replace=TRUE)
    }))
  }

  # get observations for each transect
  obs_per_sample <- by(bootdat_samps, bootdat_samps$Sample.Label,
                       function(x) unique(x$object))
  # restrict to those selected above
  obs_per_sample <- obs_per_sample[samples]

  # resample if that's what we're doing, else we keep all observations
  obs <- unlist(obs_per_sample)
  if(obs_label %in% our_resamples){
    obs <- unlist(lapply(obs_per_sample, function(x){
      sample(x, length(x), replace=TRUE)
    }))
  }

  # get all the observations (this will include NAs)
  #this_resample <- bootdat[bootdat$object %in% obs, ]
  this_resample <- bootdat[rep(1, length(obs)), ]
  for(i in 1:nrow(this_resample)){
    this_resample[i, ] <- bootdat[bootdat$object == obs[i], ]
  }

  # get the object names to be unique
  this_resample[[obs_label]] <- 1:nrow(this_resample)

  return(this_resample)
}
