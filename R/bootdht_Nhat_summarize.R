#' Simple summary for bootstrap model
#'
#' When using [`bootdht`][bootdht] one needs to use a summary function to
#' extract results from the resulting models per replicate. This function is
#' the simplest possible example of such a function, that just extracts the
#' estimated abundance (with stratum labels).
#'
#' Further examples of such functions can be found at
#' <http://examples.distancesampling.org>.
#'
#' @param ests output from [`dht2`][dht2].
#' @param fit fitted detection function object (unused).
#' @return `data.frame` with two columns ("`Nhat`" and "`Label`"), giving the
#' estimate(s) of abundance of individuals per stratum from each bootstrap
#' replicate. This `data.frame` can be examined for example, with
#' [`quantile`][stats::quantile] to compute confidence intervals.
#' @export
#' @seealso [`bootdht`][bootdht] which this function is to be used with.
bootdht_Nhat_summarize <- function(ests, fit) {
  return(data.frame(Label = ests$individuals$N$Label,
                    Nhat  = ests$individuals$N$Estimate))
}
