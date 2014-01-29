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
  cp <- cutpoints

  # pull out the distances (removing the NAs for now)
  d <- data$distance[!is.na(data$distance)]

  # check to see if any of the distances lie outside of the cutpoints
  if(any(d<cp[1]) | any(d>=cp[length(cp)])){
    stop("Some distances lie outside of the binning cutpoints. Remove and then re-try analysis.")
  }

  distbegin<-rep(NA,length(d))
  distend<-rep(NA,length(d))

  for(i in 1:(length(cp)-1)){
    # which elements of d lie between cutpoints i and i+1
    ind <- which(d>=cp[i] & d<cp[i+1])

    distbegin[ind] <- cp[i]
    distend[ind]   <- cp[i+1]
  }

  # handle NA distances, that we need to preserve
  distbegin.na <- rep(NA,length(data$distance))
  distend.na <- rep(NA,length(data$distance))
  distbegin.na[!is.na(data$distance)] <- distbegin
  distend.na[!is.na(data$distance)] <- distend

  # put all that together and make a data.frame
  data <- cbind(data,
                distbegin=distbegin.na,
                distend=distend.na)
  data <- data.frame(data)

  return(data)
}
