#' Detection function objects when detection is certain
#'
#' Create a detection function object for strip/plot surveys for use with
#' `dht2`.
#'
#' @export
#' @param data as specified for `ds` and `ddf` (including a `size` column)
#' @param width right truncation
#' @param left left truncation (default 0, no left truncation)
#' @param transect `"line"` or `"point"` transect
#' @author David L Miller
dummy_ddf <- function(data, width, left=0, transect="line"){


  if(!(transect %in% c("line", "point"))){
    stop("transect should be \"line\" or \"point\"")
  }

  df_obj <- list()

  # put object IDs in a data.frame...
  df_obj$data <- data
  object <- data$object

  # set the fitted values
  df_obj$fitted <- rep(1, length(object))
  names(df_obj$fitted) <- object

  # truncation(s)
  df_obj$meta.data <- list()
  df_obj$meta.data$width <- width
  df_obj$meta.data$left <- left

  df_obj$meta.data$point <- FALSE
  if(transect == "point"){
    df_obj$meta.data$point <- TRUE
  }

  # make the method be "dummy"
  df_obj$method <- "dummy"

  # hessian?
  df_obj$hessian <- matrix(1, 1, 1)
  df_obj$par <- 0

  class(df_obj) <- c("fake_ddf", "ds", "ddf")
  return(df_obj)
}

#' Prediction for fake detection functions
#'
#' Prediction function for dummy detection functions. The function returns as
#' many 1s as there are rows in \code{newdata}. If \code{esw=TRUE} then the
#' strip width is returned.
#'
#' @export
#' @param object model object
#' @param newdata how many 1s should we return?
#' @param compute unused, compatibility with [`mrds::predict`][mrds::predict]
#' @param int.range unused, compatability with [`mrds::predict`][mrds::predict]
#' @param esw should the strip width be returned?
#' @param \dots for S3 consistency
#' @author David L Miller
predict.fake_ddf <- function(object, newdata=NULL, compute=FALSE,
                             int.range=NULL, esw=FALSE, ...){

  ret <- list()

  if(is.null(newdata)){
    newdata <- data.frame(dummy=object$fitted)
  }

  if(esw){
    ret$fitted <- rep(object$meta.data$width-object$meta.data$left,
                      nrow(newdata))
  }else{
    ret$fitted <- rep(1, nrow(newdata))
  }

  return(ret)
}

#' @export
print.fake_ddf <- function(x, ...){
  print(summary(x))
}

#' @export
summary.fake_ddf <- function(object, ...){
  object$average.p.se <- matrix(0, 1, 1)
  class(object) <- "summary.fake_ddf"
  return(object)
}

#' @export
print.summary.fake_ddf <- function(x, ...){
  cat("\nSummary for dummy ds object \n")
  cat("Number of observations : ", nrow(x$data),"\n")
  cat("Distance range         : ", x$meta.data$left, " - ",
                                   x$meta.data$width,"\n")
  cat("\nModel : No detection function, strip transect\n\n")
  cat("AIC   : NA\n")

  invisible()
}


