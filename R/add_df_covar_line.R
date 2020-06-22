#' Add covariate levels detection function plots
#'
#' @inherit mrds::add_df_covar_line
#' @name add_df_covar_line
#' @docType methods
#' @param ddf a fitted detection function object
#' @param data a \code{data.frame} with the covariate combination you want to plot
#' @param \dots extra arguments to give to \code{\link{line}} (\code{lty}, \code{lwd}, \code{col})
#' @param ndist number of points to evaluate the detection function at
#' @note This function is located in the \code{mrds} package but the documentation is provided here for easy access.
#' @examples
#' \dontrun{
#' # example using a model for the minke data
#' data(minke)
#' # fit a model
#' result <- ds(minke, formula=~Region.Label)
#'
#' # make a base plot, showpoints=FALSE makes the plot less busy
#' plot(result, showpoints=FALSE)
#'
#' # add lines for sex one at a time
#' add_df_covar_line(result, data.frame(Region.Label="South"), lty=2)
#' add_df_covar_line(result, data.frame(Region.Label="North"), lty=3)
#'
#' # add a legend
#' legend(1.5, 1, c("Average", "South", "North"), lty=1:3)
#'
#' }
NULL
