#' Summary of distance sampling analysis 
#'
#' Provides a brief summary of a distance sampling analysis. This includes
#' parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error.
#'
#' @aliases summary.dsmodel
#' @param object a distance analysis
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return list of extracted and summarized objects
#' @note This function just calls [`summary.ds`][summary.ds] and [`dht`][dht],
#' collates and prints the results in a nice way.
#' @author David L. Miller
#' @keywords utility
#' @export
summary.dsmodel <- function(object,...){

  #  se if TRUE, computes standard errors
  #  N if TRUE, computes abundance in covered (sampled) region
  ans <- list(ds=summary(object$ddf, se=TRUE, N=TRUE),
              dht=object$dht,
              ddf=object$ddf)

  class(ans) <- "summary.dsmodel"
  return(ans)
}
