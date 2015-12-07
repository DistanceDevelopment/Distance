#' Goodness of fit testing and quantile-quantile plots
#'
#' Computes goodness of fit tests for the detection function. For binned distances this is only chi-squared. For exact distances, Kolmogorov-Smirnov and Cramer-von Mises goodness of fit tests are computed (if \code{chisq=TRUE} then chi-squared is also computed). A quantile-quantile plot is for the fitted model is produced as a graphical representation of goodness of fit (this can be suppressed by setting \code{plot=FALSE}).
#'
#' See \code{\link{ddf.gof}} for further details.
#'
#' @param model a fitted detection function.
#' @param plot if \code{TRUE} the Q-Q plot is plotted
#' @param chisq if \code{TRUE} then chi-squared statistic is calculated even for models that use exact distances. Ignored for models that use binned distances
#' @param ... other arguments to be passed to \code{\link{ddf.gof}}
#' @export
gof_ds <- function(model, plot=TRUE, chisq=FALSE, ...){

  gof <- ddf.gof(model$ddf, qq=plot, ...)

  if(model$ddf$meta.data$binned | chisq){
    return(gof)
  }else{
    gof$chisquare <- NULL
    return(gof)
  }
}
