#' Plot a fitted detection function
#'
#' This is just a simple wrapper around \code{\link{plot.ds}}.
#'
#' @param x an object of class \code{dsmodel}.
#' @return \code{NULL}, just produces a plot.
#' @S3method plot dsmodel
#' @aliases plot.dsmodel
#' @export
#' @author David L. Miller
plot.dsmodel <- function(x,...){

  plot(x$ddf,...)

  invisible()
}
