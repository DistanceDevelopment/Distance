# safely truncate a flatfile dataset
# need to ensure that:
#  1. NAs are preserved
#  2. truncation doesn't cause loss of samples
safetruncate <- function(flatfile, right, left){

  # check Sample.Labels don't get dropped
  sl <- unique(flatfile$Sample.Label)
  # indices to keep
  find <- flatfile$distance <= right &
          flatfile$distance >= left
  find[is.na(find)] <- TRUE

  fsl <- unique(flatfile$Sample.Label[find])

  # would we drop Sample.Labels?
  if(length(fsl) != length(sl)){
    # if so, get the ones we would drop
    sl_diff <- setdiff(sl, fsl)
    # add them to the keep list
    find[flatfile$Sample.Label %in% sl_diff] <- TRUE
    # set the observation-specific data to NA
    flatfile[flatfile$Sample.Label %in% sl_diff, ]$distance <- NA

    # if size or object columns are present set their values to NA
    if(!is.null(flatfile$object)){
      flatfile[flatfile$Sample.Label %in% sl_diff, ]$object <- NA
    }
    if(!is.null(flatfile$size)){
      flatfile[flatfile$Sample.Label %in% sl_diff, ]$size <- NA
    }
  }

  # keep only these rows
  flatfile <- flatfile[find, , drop=FALSE]

  return(flatfile)
}
