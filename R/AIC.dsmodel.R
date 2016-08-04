#' Akaike's An Information Criterion for detection functions
#'
#' Extract the AIC from a fitted detection function.
#'
#' @param object a fitted detection function object
#' @param k penalty per parameter to be used; the default \code{k = 2} is the "classical" AIC
#' @param \dots required for S3 but ignored
#' @author David L Miller
#' @export
#' @importFrom stats logLik
AIC.dsmodel <- function(object, ..., k=2){

  # see also logLik.dsmodel
  ll <- logLik(object)
  return(-2*ll + k*attr(ll, "df"))
}
