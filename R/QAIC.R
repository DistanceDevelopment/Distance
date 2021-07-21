#' QIAC Information Criterion for detection functions when data is overdispersed
#'
#' Overdispersion causes AIC to select overly-complex models, so analysts
#' should specify the number/order of adjustment terms manually when fitting
#' distance sampling models to data from camera traps, rather than allowing
#' automated selection using AIC. Howe et al (2019) describes two methods for
#' performing model selection of distance sampling models in the face of
#' overdispersion. In both cases comparisons should be made between models
#' which use the same key function in their detection function specification.
#'
#' The first method of Howe et al (2019) employs a two-step process. First, an
#' overdisersion factor (\deqn{\hat{c}_1}{chat1}) is computed for each key
#' function family from the most complex model in each family (derived from the
#' chi-squared goodness of fit test statistic divided by its degrees of
#' freedom). If multiple detection functions are supplied to `QAIC`, then
#' `chat` will be calculated automatically. Alternatively `chat` can be
#' supplied as an argument.
#'
#' The second method (\deqn{\hat{c}_2}{chat2}) in the notation of Howe et al)
#' uses the number of distance observations recorded per independent encounter
#' between an animal and an observer. This quantity can be calculated from the
#' raw data. In camera trap surveys of solitary animals,
#' \deqn{\hat{c}_2}{chat2} would be the mean number of distance observations
#' recorded during a single pass by an animal in front of a trap. In surveys of
#' social animals employing human observers, \deqn{\hat{c}_2}{chat2} would be
#' the mean number of detected animals per detected group, and in camera trap
#' surveys of social animals \deqn{\hat{c}_2}{chat2} the mean number of
#' distance observations recorded during an encounter between a group of
#' animals and a CT. In this case one can calculate \deqn{\hat{c}_2}{chat2} and
#' use the `chat` argument to specify it for one or more models.
#'
#' @param object a fitted detection function object
#' @param chat a value of \deqn{\hat{c}}{chat} to be used in QIAC calculation
#' @param k penalty per parameter to be used; default 2
#' @param \dots additional fitted model objects.
#' @author David L Miller, based on code from Eric Rexstad and explanation from
#' Eric Howe.
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
#'
#' # using a pre-calculated chat
#' chat <- attr(QAIC(model1, model2, model4, model6), "chat")
#' QAIC(model1, chat=chat)
#' }
QAIC <- function(object, ..., chat=NULL, k=2){

  # get the models
  models <- list(object, ...)
  models$chat <- models$k <- NULL

  # if there is only one model, there is no comparison to make
  if(length(models)<2 & is.null(chat)){
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
  npar <- unlist(lapply(models, function(x) length(x$ddf$par)))
  aic <-  unlist(lapply(models, function(x) x$ddf$criterion))

  if(is.null(chat)){
    chat <- chat(models[[which.max(npar)]])
  }
  qaics <- data.frame(df=npar,
                      QAIC=unlist(lapply(models, qaic, chat=chat, k=k)))

  # add row names
  call <- match.call(expand.dots=TRUE)
  call$chat <- call$k <- NULL
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
qaic <- function(modobj, chat, k) {
  2* modobj$ddf$ds$value/chat + k * (length(modobj$ddf$pars)+1)
}
