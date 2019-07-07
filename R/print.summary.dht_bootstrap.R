#' @export
print.summary.dht_bootstrap <- function(x, ...){
object <- x
  cat("Bootstrap results\n\n")
  cat("Boostraps          :", object$nboot, "\n")
  cat("Successes          :", object$nbootsuccess, "\n")
  cat("Failures           :", object$nbootfail, "\n")

  cat("\n")

  object$tab <- round(object$tab, 2)

  print(object$tab)
}
