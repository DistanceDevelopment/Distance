#' Create bins from a set of binned distances and a set of cutpoints.
#'
#' This is a service routine and shouldn't be necessary in normal analyses.
#'
#' @param data \code{data.frame} with at least the column \code{distance}.
#' @param cutpoints vector of cutpoints for the bins
#'
#' @return data \code{data} with two extra columns \code{distbegin} and
#'        \code{distend}.
#'
#' @author David L. Miller
#'
#'
create.bins <- function(data,cutpoints){

  # lazy typist
  cp<-cutpoints

  # pull out the distances
  d<-data$distance

  # check to see if any of the distances lie outside of the cutpoints
  if(any(d<cp[1]) | any(d>cp[length(cp)])){
    stop("Some distances lie outside of the binning cutpoints. Remove and then re-try analysis.")
  }

  distbegin<-rep(NA,length(d))
  distend<-rep(NA,length(d))

  # catch anything in the first bin that was recorded as the
  # bottom cutpoint
  ind0 <- which(d==cp[1])
  distbegin[ind0] <- cp[1]
  distend[ind0] <- cp[2]

  for(i in 1:(length(cp)-1)){
    # which elements of d lie between cutpoints i and i+1
    ind <- which(d>=cp[i] & d<cp[i+1])

    distbegin[ind] <- cp[i]
    distend[ind]   <- cp[i+1]
  }

  data <- cbind(data,distbegin=distbegin,distend=distend)
  data <- data.frame(data)

  return(data)
}
