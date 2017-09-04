#' log-likelihood value for a fitted detection function
#'
#' Extract the log-likelihood from a fitted detection function.
#'
#' @param object a fitted detection function model object
#' @param \dots included for S3 completeness, but ignored
#' @return a numeric value giving the log-likelihood with two attributes: \code{"df"} the "degrees of freedom" for the model (number of parameters) and \code{"nobs"} the number of observations used to fit the model
#' @export
#' @author David L Miller
#' @examples
#' \dontrun{
#' library(Distance)
#' data(minke)
#' model <- ds(minke, truncation=4)
#' # extract the log likelihood
#' logLik(model)
#' }
logLik.dsmodel <- function(object, ...){

  # see ?logLik for information on why

  ret <- object$ddf$lnl

  attr(ret, "df") <- length(object$ddf$par)
  attr(ret, "nobs") <- nrow(object$ddf$data)

  class(ret) <- "logLik"
  return(ret)
}
