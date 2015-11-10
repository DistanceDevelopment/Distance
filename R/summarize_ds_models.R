#' Make a table of summary statistics for detection function models
#'
#' Provide a summary table of useful information about fitted detection functions. This can be useful when paired with \code{knitr}s \code{\link{kable}} function.
#'
#' Note that the column names are in LaTeX format, so if you plan to manipulate the resulting \code{data.frame} in R, you may wish to rename the columns for ease of access.
#'
#' @param ... models to be summarised
#' @param sort column to sort by (default \code{"AIC"})
#' @author David L Miller
#' @export
summarize_ds_models <- function(..., sort="AIC"){

  # get the models
  models <- list(...)

  # this function extracts the model data for a single model (row)
  extract_model_data <- function(model){
    summ <- summary(model)
    ret <- c(model.description(model$ddf),
             model$ddf$ds$aux$ddfobj$scale$formula,
             ddf.gof(model$ddf, qq=FALSE)$dsgof$CvM$p,
             summ$ds$average.p,
             summ$ds$average.p.se,
             model$ddf$criterion
            )
    return(ret)
  }

  # applying that to all the models then putting it into a data.frame
  res <- as.data.frame(t(as.data.frame(lapply(models, extract_model_data))),
                        stringsAsFactors=FALSE)

  # making sure the correct columns are numeric
  res[,3:6] <- apply(res[,3:6], 2, as.numeric)

  # giving the columns names
  colnames(res) <- c("Key function",
                     "Formula",
                     "C-vM $p$-value",
                     "$\\hat{P_a}$",
                     "se($\\hat{P_a}$)",
                     "AIC")
  # remove row names
  rownames(res) <- NULL

  # creating a new column for the AIC difference to the best model
  res[["$\\Delta$AIC"]] <- res$AIC - min(res$AIC, na.rm=TRUE)
  # ordering the model by AIC score
  res <- res[order(res[[sort]]),]

  # returning the data.frame
  return(res)
}
