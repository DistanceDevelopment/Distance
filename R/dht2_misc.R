# misc internal functions for dht2/bootdht
check_sample_fraction <- function(sample_fraction){
  if(length(sample_fraction)>1){
    stop("sample_fraction must be a single number")
  }
  if(sample_fraction <= 0){
    stop("sample_fraction must be > 0")
  }
  invisible(TRUE)
}
