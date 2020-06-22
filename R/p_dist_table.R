#' Distribution of probabilities of detection
#'
#' @name p_dist_table
#' @inherit mrds::p_dist_table
#' @param object fitted detection function
#' @param bins how the results should be binned
#' @param proportion should proportions be returned as well as counts?
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

