#' @export
#' @importFrom stats median quantile
summary.dht_bootstrap <- function(object, alpha=0.05, ...){

  x <- list()

  x$nboot <- attr(object, "nboot")
  x$nbootsuccess <- attr(object, "nboot")-attr(object, "nbootfail")
  x$nbootfail <- attr(object, "nbootfail")
  x$alpha <- alpha

  class(object) <- "list"
  object <- as.data.frame(object)

  numcols <- unlist(lapply(object, is.numeric))

  # build a summary object
  sumfun <- function(x){
    if(is.numeric(x)){
      xx <- data.frame(Median   = median(x, na.rm=TRUE),
                       Mean     = mean(x, na.rm=TRUE),
                       se       = sqrt(var(x, na.rm=TRUE)),
                       lcl      = quantile(x, (alpha/2)),
                       ucl      = quantile(x, 1-(alpha/2)))
      xx$cv <- xx$se/xx$Median
    }else{
      xx <- NULL#data.frame(Estimate = NA,
                #       se       = NA,
                #       ucl      = NA,
                #       lcl      = NA,
                #       cv       = NA)
    }
    return(xx)
  }

  tn <- lapply(object, sumfun)
  x$tab <- do.call(rbind.data.frame, tn)

  class(x) <- "summary.dht_bootstrap"
  return(x)
}
