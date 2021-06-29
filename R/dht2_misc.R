# misc internal functions for dht2/bootdht
# calculate covered area
area_calc <- function(width, effort, transect_type, sample_fraction){
  # Note - this function does not need to account for left truncation, as
  # it is taken care of when calculating average detection prob (by setting
  # detection prob to 0 between 0 and left truncation distance)
  if(transect_type=="point"){
    return(effort*pi*width^2*sample_fraction)
  }else{
    return(effort*2*width*sample_fraction)
  }
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
