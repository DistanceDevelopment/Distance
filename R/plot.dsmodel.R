#' Plot a fitted detection function
#'
#' This is just a simple wrapper around [`plot.ds`][mrds::plot.ds]. See the
#' manual page for that function for more information.
#'
#' @param x an object of class `dsmodel`.
#' @param pl.den shading density for histogram (default `0`, no shading)
#' @param ... extra arguments to be passed to [`plot.ds`][mrds::plot.ds].
#' @return `NULL`, just produces a plot.
#' @aliases plot.dsmodel
#' @export
#' @author David L. Miller
#' @importFrom graphics plot
#' @seealso [`add_df_covar_line`][add_df_covar_line]
plot.dsmodel <- function(x, pl.den=0, ...){

  plot(x$ddf, pl.den=pl.den, ...)

  invisible()
}
