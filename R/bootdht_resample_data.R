bootdht_resample_data <- function(bootdat, our_resamples,
                                  stratum_label="Region.Label",
                                  sample_label="Sample.Label",
                                  obs_label="object"){

  # unique combinations
  bf <- unique(bootdat[, c(stratum_label, sample_label, obs_label)])

  # stratum resample
  # which strata?
  strata <- unique(bootdat[[stratum_label]])
  # resample if that's what we're doing, else we keep all strata
  if(stratum_label %in% our_resamples){
    strata <- sample(strata, length(strata), replace=TRUE)
  }
  strata_df <- cbind(as.character(strata), make.unique(as.character(strata)))
  colnames(strata_df) <- c(stratum_label, ".new_stratum_label")

  # get all samples per stratum
  samps_per_strata <- unique(bf[, c(stratum_label, sample_label)])
  samps_per_strata <- by(bf[,c(stratum_label, sample_label)],
                         bf[[stratum_label]],
                         function(x) unique(x[[sample_label]]))


  # resample if that's what we're doing, else we keep all samples
  samples <- samps_per_strata
  if(sample_label %in% our_resamples){
    samples <- lapply(samps_per_strata, function(x){
      sample(x, length(x), replace=TRUE)
    })
  }

  # get observations for each transect
  obs_per_sample <- by(bf[,c(sample_label, obs_label)],
                       bf[[sample_label]],
                       function(x) unique(x[[obs_label]]))
  # restrict to those selected above
  obs_per_sample <- obs_per_sample[unlist(samples)]

  # resample if that's what we're doing, else we keep all observations
  obs <- obs_per_sample
  if(obs_label %in% our_resamples){
    obs <- lapply(obs_per_sample, function(x){
      sample(x, length(x), replace=TRUE)
    })
  }

  # remap sample labels
  remap <- make.unique(names(obs))

  # now populate each list element with the corresponding observations
  # or empty transects
  for(ii in 1:length(obs)){
    if(all(is.na(obs[[ii]]))){
      obs[[ii]] <- bootdat[bootdat[[sample_label]] == names(obs)[ii], ]
    }else{
      obs[[ii]] <- bootdat[bootdat$object %in% obs[[ii]], ]
    }
    obs[[ii]][[sample_label]] <- remap[ii]
  }

  # concatenate list elements to data.frame
  rr <- do.call(rbind.data.frame, obs)

  # reset the object IDs to be unique (where there are observations)
  rr[[obs_label]][is.na(rr[["distance"]])] <- NA
  rr[[obs_label]][!is.na(rr[["distance"]])] <- 1:length(rr[[obs_label]][!is.na(rr[["distance"]])])

  return(rr)
}
