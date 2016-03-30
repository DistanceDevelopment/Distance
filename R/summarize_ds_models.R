#' Make a table of summary statistics for detection function models
#'
#' Provide a summary table of useful information about fitted detection functions. This can be useful when paired with \code{knitr}s \code{kable} function.
#'
#' Note that the column names are in LaTeX format, so if you plan to manipulate the resulting \code{data.frame} in R, you may wish to rename the columns for ease of access.
#'
#' @param ... models to be summarised
#' @param sort column to sort by (default \code{"AIC"})
#' @param output should the output be given in \code{"latex"} compatible format or as \code{"plain"} text?
#' @param delta_only only output AIC differences (default \code{TRUE})
#' @author David L Miller
#' @export
summarize_ds_models <- function(..., sort="AIC", output="latex", delta_only=TRUE){

  # get the models
  models <- list(...)
  # get the model names
  model_names <- setdiff(as.character(match.call(expand.dots=TRUE)),
                         as.character(match.call(expand.dots=FALSE)))

  # this function extracts the model data for a single model (row)
  extract_model_data <- function(model){
    summ <- summary(model)

    # handle (uniform) no formula case
    formula <- model$ddf$ds$aux$ddfobj$scale$formula
    if(is.null(formula)) formula <- NA

    desc <- gsub(" key function","",model.description(model$ddf))
    ret <- c(desc,
             formula,
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

  if(output == "latex"){
    model_names <- gsub("_", '\\\\char`_', model_names)
    model_names <- paste0("\\texttt{", model_names, "}")
    res <- cbind.data.frame(model_names, res)
  }else if(output=="plain"){
    res <- cbind.data.frame(model_names, res)
  }else{
    stop("Invalid output format")
  }

  # making sure the correct columns are numeric
  res[,4:7] <- apply(res[,4:7], 2, as.numeric)

  # giving the columns names
  if(output == "latex"){
    colnames(res) <- c("Model",
                       "Key function",
                       "Formula",
                       "C-vM $p$-value",
                       "$\\hat{P_a}$",
                       "se($\\hat{P_a}$)",
                       "AIC")
  }else if(output=="plain"){
    colnames(res) <- c("Model",
                       "Key function",
                       "Formula",
                       "C-vM p-value",
                       "Average detectability",
                       "se(Average detectability)",
                       "AIC")
  }else{
    stop("Invalid output format")
  }
  # remove row names
  rownames(res) <- NULL

  # creating a new column for the AIC difference to the best model
  if(output == "latex"){
    res[["$\\Delta$AIC"]] <- res$AIC - min(res$AIC, na.rm=TRUE)
  }else if(output=="plain"){
    res[["Delta AIC"]] <- res$AIC - min(res$AIC, na.rm=TRUE)
  }
  # ordering the model by AIC score
  res <- res[order(res[[sort]]),]

  # remove the AIC column if asked and just return the deltas
  if(delta_only){
    res$AIC <- NULL
  }

  # returning the data.frame
  return(res)
}
