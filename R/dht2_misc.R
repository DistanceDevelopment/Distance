# misc internal functions for dht2/bootdht
# calculate covered area
area_calc <- function(width, left, effort, transect_type, sample_fraction){

  res <- rep(NA, length(width))

  res[transect_type=="point"] <- effort[transect_type=="point"]*pi*
                                  width[transect_type=="point"]^2*
                                  sample_fraction[transect_type=="point"]
  res[transect_type=="line"] <- effort[transect_type=="line"]*2*
                                 width[transect_type=="line"]*
                                 sample_fraction[transect_type=="line"]

  return(res)
}

# sample fractions
dht2_sample_fraction <- function(sample_fraction, data){

  # do some checking
  if(!is.data.frame(sample_fraction) && length(sample_fraction)>1){
    stop("sample_fraction must be a single number or a data.frame")
  }
  if(!is.data.frame(sample_fraction) && sample_fraction <= 0){
    stop("sample_fraction must be > 0")
  }

  if(is.data.frame(sample_fraction) &&
     !all(names(sample_fraction) %in% c("Sample.Label", "fraction"))){
    stop("sample_fraction data.frame columns must be \"Sample.Label\" and \"fraction\"")
  }

  if(is.data.frame(sample_fraction)){
    sample_fraction$sample_fraction <- sample_fraction$fraction
    sample_fraction$fraction <- NULL
    data <- merge(data, sample_fraction, by="Sample.Label")
  }else{
    data$sample_fraction <- sample_fraction
  }

  return(data)
}
