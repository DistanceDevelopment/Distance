#' Create bins from a set of binned distances and a set of cutpoints.
#'
#' `create.bins` is now deprecated, please use [`create_bins`][create_bins]
#'
#' @param data `data.frame` with at least the column `distance`.
#' @param cutpoints vector of cutpoints for the bins
#'
#' @return argument `data` with two extra columns `distbegin` and
#'        `distend`.
#'
#' @author David L. Miller
#' @export
create.bins <- function(data, cutpoints){
  stop("create.bins is deprecated, please use create_bins")
}
