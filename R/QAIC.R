#' Tools for model selection when distance sampling data are overdispersed
#'
#' Overdispersion causes AIC to select overly-complex models, so analysts
#' should specify the number/order of adjustment terms manually when fitting
#' distance sampling models to data from camera traps, rather than allowing
#' automated selection using AIC. Howe et al (2019) described a two-step method
#' for selecting among models of the detection function in the face of
#' overdispersion.
#'
#' In step 1, and overdispersion factor (\eqn{\hat{c}}{chat}) is computed
#' either (1) for each key function family, from the most complex model in each
#' family, as the chi-square goodness of fit test statistic divided by its
#' degrees of freedom (\eqn{\hat{c}_1}{chat1}), or (2) for all models in the
#' candidate set, from the raw data (\eqn{\hat{c}_1}{chat2}). In camera trap
#' surveys of solitary animals, \eqn{\hat{c}_1}{chat2} would be the mean number
#' of distance observations recorded during a single pass by an animal in front
#' of a trap. In surveys of social animals employing human observers,
#' \eqn{\hat{c}_1}{chat2} would be the mean number of detected animals per
#' detected group, and in camera trap surveys of social animals
#' \eqn{\hat{c}_1}{chat2} the mean number of distance observations recorded
#' during an encounter between a group of animals and a CT.  In step two, the
#' chi-square goodness of fit statistic divided by its degrees of freedom is
#' calculated for the QAIC-minimizing model within each key function, and the
#' model with the lowest value is selected for estimation.
#'
#' The `QAIC()` function should only be used select among models with the same
#' key function (step 1). `QAIC()` uses \eqn{\hat{c}_1}{chat1} by default,
#' computing it from the model with the most parameters. Alternatively,
#' \eqn{\hat{c}_1}{chat2} can be calculated from the raw data and included in
#' the call to `QAIC()`. Users must identify the QAIC-minimizing model within
#' key functions from the resulting `data.frame`, and provide these models as
#' arguments in `ch2_select()`. `chi2_select()` then computes and reports the
#' chi-square goodness of fit statistic divided by its degrees of freedom for
#' each of those models. The model with the lowest value is recommended for
#' estimation.
#'
#' @param object a fitted detection function object
#' @param chat a value of \eqn{\hat{c}}{chat} to be used in QAIC calculation
#' @param k penalty per parameter to be used; default 2
#' @param \dots additional fitted model objects.
#' @return a `data.frame` with one row per model supplied, in the same order as
#'         given
#' @author David L Miller, based on code from Eric Rexstad and explanation from
#' Eric Howe.
#' @references Howe, E. J., Buckland, S. T., Després-Einspenner, M.-L., & Kühl, H. S. (2019). Model selection with overdispersed distance sampling data. Methods in Ecology and Evolution, 10(1), 38–47. \doi{10.1111/2041-210X.13082}
#' @export
#' @importFrom stats logLik
#' @aliases chi2_select
#' @examples
#' \dontrun{
#' library(Distance)
#' data("wren_cuecount")
#'
#' # fit hazard-rate key models
#' w3.hr0 <- ds(wren_cuecount, transect="point", key="hr", adjustment=NULL,
#'              truncation=92.5)
#' w3.hr1 <- ds(wren_cuecount, transect="point", key="hr", adjustment="cos",
#'              order=2, truncation=92.5)
#' w3.hr2 <- ds(wren_cuecount, transect="point", key="hr", adjustment="cos",
#'              order=c(2, 4), truncation=92.5)
#'
#' # fit unform key models
#' w3.u1 <- ds(wren_cuecount, transect="point", key="unif", adjustment="cos",
#'             order=1, truncation=92.5)
#' w3.u2 <- ds(wren_cuecount, transect="point", key="unif", adjustment="cos",
#'             order=c(1,2), monotonicity="none",  truncation=92.5)
#' w3.u3 <- ds(wren_cuecount, transect="point", key="unif", adjustment="cos",
#'             order=c(1,2,3), monotonicity="none", truncation=92.5)
#'
#' # fit half-normal key functions
#' w3.hn0 <- ds(wren_cuecount, transect="point", key="hn", adjustment=NULL,
#'              truncation=92.5)
#' w3.hn1 <- ds(wren_cuecount, transect="point", key="hn", adjustment="herm",
#'              order=2, truncation=92.5)
#'
#' # stage 1: calculate QAIC per model set
#' QAIC(w3.hr0, w3.hr1, w3.hr2)  # no adjustments smallest
#' QAIC(w3.u1, w3.u2, w3.u3)     # 2 adjustment terms (by 0.07)
#' QAIC(w3.hn0, w3.hn1)  # no adjustments smallest
#'
#' # stage 2: select using chi^2/degrees of freedom between sets
#' chi2_select(w3.hr0, w3.u2, w3.hn0)
#'
#' # example using a pre-calculated chat
#' chat <- attr(QAIC(w3.hr0, w3.hr1, w3.hr2), "chat")
#' QAIC(w3.hr0, chat=chat)
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
  # c-hat is computed for the most parameter-rich model in the group QAIC is
  # calculated for each model in group based upon this c-hat

  # check all models have the same key function
  keys <- unlist(lapply(models, function(x) x$ddf$ds$aux$ddfobj$type))
  if(length(unique(keys))!=1){
    stop("All key functions must be the same")
  }

  num.models <- length(models)
  npar <- unlist(lapply(models, function(x) length(x$ddf$par)))

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

#' @export
#' @rdname QAIC
chi2_select <- function(object, ...){

  # get the models
  models <- list(object, ...)

  # if there is only one model, there is no comparison to make
  if(length(models)<2 & is.null(chat)){
    stop("Only 1 model specified, no model selection can be performed")
  }

  # check all models have the same key function
  keys <- unlist(lapply(models, function(x) x$ddf$ds$aux$ddfobj$type))
  if(length(unique(keys))!= length(keys)){
    stop("All key functions must be different")
  }

  num.models <- length(models)


  # add chi^2/df
  res <- data.frame(criteria = unlist(lapply(models, chat)))
  # add row names
  call <- match.call(expand.dots=TRUE)
  rownames(res) <- as.character(call)[-1]
  # add chat
  attr(res, "chat") <- chat

  res
}

# compute c-hat for a dsmodel object using Method 1 of Howe et al. (2018)
chat <- function(modobj) {
  test <- gof_ds(modobj, plot=FALSE, chisq=TRUE)
  test$chisquare$chi1$chisq / test$chisquare$chi1$df
}

# compute QAIC for a dsmodel object given a c-hat
qaic <- function(modobj, chat, k) {
  - 2* modobj$ddf$lnl/chat + k * (length(modobj$ddf$par)+1)
}
