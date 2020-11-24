#' Add covariate levels detection function plots
#'
#' @inherit mrds::add_df_covar_line
#' @name add_df_covar_line
#' @docType methods
#' @param ddf a fitted detection function object.
#' @param data a \code{data.frame} with the covariate combination you want to plot.
#' @param \dots extra arguments to give to \code{\link{line}} (\code{lty}, \code{lwd}, \code{col}).
#' @param ndist number of distances at which to evaluate the detection function.
#' @param pdf should the line be drawn on the probability density scale; ignored for line transects
#' @param breaks required to ensure that PDF lines are the right size, should match what is supplied to original \code{plot} command. Defaults to "Sturges" breaks, as in \code{\link{hist}}. Only used if \code{pdf=TRUE}
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
#' # point transect example
#' data(amakihi)
#' result <- ds(amakihi, truncation=150, transect="point", formula=~OBs)
#' plot(result, showpoints=FALSE, pdf=TRUE)
#' add_df_covar_line(result,
#'                   data.frame(OBs=na.omit(unique(amakihi$OBs))), pdf=TRUE)
#' }
NULL
