#' Simple summary of density results for bootstrap model
#'
#' When using [`bootdht`][bootdht] one needs to use a summary function to
#' extract results from the resulting models per replicate. This function is
#' the simplest possible example of such a function, that just extracts the
#' estimated density (with stratum labels).
#'
#' Further examples of such functions can be found at
#' <http://examples.distancesampling.org>.
#'
#' @param ests output from [`dht2`][dht2].
#' @param fit fitted detection function object (unused).
#' @return `data.frame` with two columns ("`Dhat`" and "`Label`"), giving the
#' estimate(s) of density of individuals per stratum from each bootstrap
#' replicate. This `data.frame` can be examined for example, with
#' [`quantile`][stats::quantile] to compute confidence intervals.
#' @export
#' @seealso [`bootdht`][bootdht] which this function is to be used with and
#' [`bootdht_Nhat_summarize`][bootdht_Nhat_summarize] which does the same job
#' but returns abundance results.
bootdht_Dhat_summarize <- function(ests, fit) {
  return(data.frame(Label = ests$individuals$D$Label,
                    Dhat  = ests$individuals$D$Estimate))
}
