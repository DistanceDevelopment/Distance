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
#' @examples
#' \dontrun{
#' # fit and test a simple model for the golf tee data
#' library(Distance)
#' data(book.tee.data)
#' tee.data<-book.tee.data$book.tee.dataframe[book.tee.data$book.tee.dataframe$observer==1,]
#' ds.model <- ds(tee.data,4)
#' # don't make plot
#' gof_ds(ds.model, plot=FALSE)
#'}
gof_ds <- function(model, plot=TRUE, chisq=FALSE, ...){

  gof <- suppressMessages(ddf.gof(model$ddf, qq=plot, ...))

  if(model$ddf$meta.data$binned | chisq){
    return(gof)
  }else{
    gof$chisquare <- NULL
    return(gof)
  }
}
