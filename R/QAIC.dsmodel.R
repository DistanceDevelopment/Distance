#' QIAC Information Criterion for detection functions when data is overdispersed
#'
#' Overdispersion causes AIC to select overly-complex models, so analysts
#' should specify the number/order of adjustment terms manually when fitting
#' distance sampling models to data from camera traps, rather than allowing
#' automated selection using AIC. Howe et al (2019) describes two methods for
#' performing model selection of distance sampling models in the face of
#' overdispersion.
#'
#' The first method of Howe et al (2019) employs a two-step process. First, an
#' overdisersion factor (`chat`) is computed for each key function family from
#' the most complex model in each family (derived from the
#' chi-squared goodness of fit test statistic divided by its degrees of
#' freedom).
#'
#' @param object a fitted detection function object
#' @param \dots optionally more fitted model objects.
#' @author David L Miller, based on code from Eric Rexstad
#' @references Howe, E. J., Buckland, S. T., Després-Einspenner, M.-L., & Kühl, H. S. (2019). Model selection with overdispersed distance sampling data. Methods in Ecology and Evolution, 10(1), 38–47. \doi{10.1111/2041-210X.13082}
#' @export
#' @importFrom stats logLik
#' @examples
#' \dontrun{
#' library(Distance)
#' data(minke)
#' model1 <- ds(minke, truncation=4, adjustment=NULL)
#' model2 <- ds(minke, truncation=4, adjustment="cos", order=2)
#' model4 <- ds(minke, truncation=4, adjustment="cos", order=c(2, 4))
#' model6 <- ds(minke, truncation=4, adjustment="cos", order=c(2, 4, 6))
#' # calculate QAIC
#' QAIC(model1, model2, model4, model6)
# }
QAIC.dsmodel <- function(object, ...){

  # get the models
  models <- list(object, ...)

  # if there is only one model, there is no comparison to make
  if(length(models)<2){
    stop("Only 1 model specified, no model selection can be performed")
  }

  # based on qaic.pass1
  # Performs Pass 1 model selection based upon Method 1 of Howe et al. (2018)
  # c-hat is computed for the most parameter-rich model in the group qaic is
  # calculated for each model in group based upon this c-hat

  # check all models have the same key function
  keys <- unlist(lapply(models, function(x) x$ddf$ds$aux$ddfobj$type))
  if(length(unique(keys))!=1){
    stop("All key functions must be the same")
  }

  num.models <- length(models)
  npar <- unlist(lapply(models, function(x) length(x$ddf$ds$par)))
  aic <-  unlist(lapply(models, function(x) x$ddf$criterion))

  chat.bigmod <- chat(models[[which.max(npar)]])
  qaics <- data.frame(df=npar,
                      QAIC=unlist(lapply(models, qaic, chat=chat.bigmod)))

  # add row names
  call <- match.call(expand.dots=TRUE)
  rownames(qaics) <- as.character(call)[-1]

  # add chat
  attr(qaics, "chat") <- chat

  qaics
}

# compute c-hat for a dsmodel object using Method 1 of Howe et al. (2018)
chat <- function(modobj) {
  test <- gof_ds(modobj, plot=FALSE, chisq=TRUE)
  test$chisquare$chi1$chisq / test$chisquare$chi1$df
}

# compute QAIC for a dsmodel object given a c-hat
qaic <- function(modobj, chat) {
  2* modobj$ddf$ds$value/chat + 2 * (length(modobj$ddf$ds$pars)+1)
}
