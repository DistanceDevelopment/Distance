#' Add covariate levels detection function plots
#'
#' @inherit mrds::add_df_covar_line parameters return references description details sections seealso
#' @name add_df_covar_line
#' @docType methods
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
