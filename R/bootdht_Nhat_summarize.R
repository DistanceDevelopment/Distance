#' Simple summary for bootstrap model
#'
#' When using \code{\link{bootdht}} one needs to use a summary function to extract results from the resulting models per replicate. This function is the simplest possible example of such a function, that just extracts the estimated abundance.
#' Further examples of such functions can be found at \url{http://examples.distancesampling.org}.
#'
#' @param ests output from \code{\link{dht2}}.
#' @param fit fitted detection function object (unused).
#' @return \code{data.frame} with one column ("\code{Nhat}"), containing estimate(s) of abundance of individuals from each bootstrap replicate.  This data frame can be examined for example, with \code{quantile} to compute confidence intervals.
#' @export
#' @seealso \code{\link{bootdht}}, which this function is to be used with.
bootdht_Nhat_summarize <- function(ests, fit) {
  return(data.frame(Nhat=ests$individuals$N$Estimate))
}
