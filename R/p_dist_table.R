#' Distribution of probabilities of detection
#'
#' @inherit mrds::add_df_covar_line param return references description details sections seealso
#' @name p_dist_table
#' @docType methods
#' @note This function is located in the \code{mrds} package but the documentation is provided here for easy access.
#' @examples
#' \dontrun{
#' # example using a model for the minke data
#' data(minke)
#' # fit a model
#' result <- ds(minke, formula=~Region.Label)
#' # print table
#' p_dist_table(result)
#' # with proportions
#' p_dist_table(result, proportion=TRUE)
#' }
NULL

