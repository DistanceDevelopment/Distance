#' Make a table of summary statistics for detection function models
#'
#' Provide a summary table of useful information about fitted detection functions. This can be useful when paired with \code{knitr}s \code{kable} function. By default models are sorted by AIC and will therefore not allow models with different truncations and distance binning.
#'
#' Note that the column names are in LaTeX format, so if you plan to manipulate the resulting \code{data.frame} in R, you may wish to rename the columns for ease of access.
#'
#' @param ... models to be summarised
#' @param sort column to sort by (default \code{"AIC"})
#' @param output should the output be given in \code{"latex"} compatible format or as \code{"plain"} text?
#' @param delta_only only output AIC differences (default \code{TRUE})
#' @author David L Miller
#' @export
#' @examples
#' \dontrun{
#' # fit some models to the golf tee data
#' library(Distance)
#' data(book.tee.data)
#' tee.data<-book.tee.data$book.tee.dataframe[book.tee.data$book.tee.dataframe$observer==1,]
#' model_hn <- ds(tee.data,4)
#' model_hr <- ds(tee.data,4, key="hr")
#' summarize_ds_models(model_hr, model_hn, output="plain")
#'}
summarize_ds_models <- function(..., sort="AIC", output="latex", delta_only=TRUE){

  # get the models
  models <- list(...)

  # get the model names
  model_names <- setdiff(as.character(match.call(expand.dots=TRUE)),
                         as.character(match.call(expand.dots=FALSE)))


  ## checking
  # can't compare models with different truncations
  r_truncs <- unlist(lapply(models, function(x) x$ddf$meta.data$width))
  l_truncs <- unlist(lapply(models, function(x) x$ddf$meta.data$left))
  if(!all(abs(c(r_truncs-mean(r_truncs),
                l_truncs-mean(l_truncs))) < 1e-8)){
    stop("All truncations must be the same for AIC comparison.")
  }
  # check all binned
  binned <- unlist(lapply(models, function(x) x$ddf$meta.data$binned))
  if((any(binned) & !all(binned)) | (any(!binned) & !all(!binned))){
    stop("Can't compare binned and unbinned distances")
  }
  # check all binning is the same
  if(all(binned)){
    breaks <- lapply(models, function(x) x$ddf$meta.data$breaks)
    # if the breaks aren't the same length it's easy
    len_breaks <- unlist(lapply(breaks, length))
    if(!all(abs(len_breaks-mean(len_breaks)) < 1e-8)){
      stop("Distance binning must be the same for all models.")
    }
    # if not??? (WARNING: Byzantine process :( )
    for(i in seq_along(breaks[[1]])){
      this_set <- unlist(lapply(breaks, "[[", i))
      if(!all(abs(this_set-mean(this_set)) < 1e-8)){
        stop("Distance binning must be the same for all models.")
      }
    }
  }

  # this function extracts the model data for a single model (row)
  extract_model_data <- function(model){
    summ <- summary(model)

    # handle (uniform) no formula case
    formula <- model$ddf$ds$aux$ddfobj$scale$formula
    if(is.null(formula)) formula <- NA

    desc <- gsub(" key function","",model.description(model$ddf))
    # only get CvM if not binned
    if(model$ddf$meta.data$binned){
      gof <- suppressMessages(gof_ds(model, chisq=TRUE)$chisquare$chi1$p)
    }else{
      gof <- suppressMessages(ddf.gof(model$ddf, qq=FALSE)$dsgof$CvM$p)
    }
    ret <- c(desc,
             formula,
             gof,
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

  # what test did we do?
  if(all(binned)){
    gof_name <- "Chi^2 p-value"
    gof_latexname <- "$\\chi^2$ $p$-value"
  }else{
    gof_name <- "C-vM $p$-value"
    gof_latexname <- "C-vM p-value"
  }

  # giving the columns names
  if(output == "latex"){
    colnames(res) <- c("Model",
                       "Key function",
                       "Formula",
                       gof_latexname,
                       "$\\hat{P_a}$",
                       "se($\\hat{P_a}$)",
                       "AIC")
  }else if(output=="plain"){
    colnames(res) <- c("Model",
                       "Key function",
                       "Formula",
                       gof_name,
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
