#' Akaike's An Information Criterion for detection functions
#'
#' Extract the AIC from a fitted detection function.
#'
#' @param object a fitted detection function object
#' @param k penalty per parameter to be used; the default \code{k = 2} is the "classical" AIC
#' @param \dots optionally more fitted model objects.
#' @author David L Miller
#' @export
#' @importFrom stats logLik
#' @examples
#' \dontrun{
#' library(Distance)
#' data(minke)
#' model <- ds(minke, truncation=4)
#' model_hr <- ds(minke, truncation=4, key="hr")
#' # extract the AIC for 2 models
#' AIC(model, model_hr)
#' }
AIC.dsmodel <- function(object, ..., k=2){

  # get the models
  models <- list(object, ...)
  models$k <- NULL

  # build the table
  aics <- matrix(NA, nrow=length(models), ncol=2)
  for(i in seq_along(models)){
    ll <- logLik(models[[i]])
    aics[i, 1] <- attr(ll, "df")
    aics[i, 2] <- -2*ll + k*attr(ll, "df")
  }
  # make it a data.frame
  aics <- as.data.frame(aics)
  names(aics) <- c("df", "AIC")
  # add row names
  call <- match.call(expand.dots=TRUE)
  call$k <- NULL
  rownames(aics) <- as.character(call)[-1]

  return(aics)
}
