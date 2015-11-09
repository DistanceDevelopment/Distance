#' Goodness of fit testing and quantile-quantile plots
#'
#' Computes chi-squared, Kolmogorov-Smirnov and Cramer-von Mises goodness of fit tests for the detection function. Optionally plot a quantile-quantile plot for fitted model as a graphical representation of goodness of fit.
#'
#' See \code{\link{ddf.gof}} for details.
#'
#' @param model a fitted detection function.
#' @param plot if \code{TRUE} the Q-Q plot is plotted
#' @export
gof_ds <- function(model, plot=TRUE, ...){
  ddf.gof(model$ddf, qq=plot, ...)
}
