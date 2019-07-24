# internal function to check input data
dht2_checkdata <- function(ddf, observations, transects, geo_strat, strat_formula,
                           stratum_labels, geo_stratum_labels){

  # required columns for observations and transects
  col_req <- list(observations = c("object", "Sample.Label"),
                  transects    = c("Sample.Label", "Effort"))
  # if we have geographical stratification, add Area too
  if(!is.null(geo_strat)){
    col_req[["geo_strat"]] <- c(geo_stratum_labels, "Area")
  }

  # check columns are right
  # loop over the data frames
  for(dname in names(col_req)){
    # first check for required columns
    if(!all(col_req[[dname]] %in% colnames(get(dname)))){
      stop(paste0(dname, " data must at least contain columns named ",
                  paste0(col_req[[dname]], collapse=", ")))
    }

    # then check for the additional stratum names
    # get the names of all columns
    these_colnames <- unlist(lapply(names(col_req), function(x) colnames(get(x))))

    # check that the stratum labels are *somewhere*
    # if they are remove from this list
    stratum_labels <- stratum_labels[!(stratum_labels %in% these_colnames)]
  }

  # check the stratum label isn't in the ddf data
  stratum_labels <- stratum_labels[!(stratum_labels %in% colnames(ddf$data))]

  # are there any stratum labels left? If so, something bad has happened
  if(length(stratum_labels) > 0){
    stop(paste0("Stratification variable(s) \"",
                stratum_labels,
                "\" not found in the data"))
  }

  # check that Areas and labels are consistent over the geographical strata
  if(!identical(geo_strat, unique(geo_strat))){
    stop("Inconsistent Areas/stratum labels in `geo_strat`")
  }
}
