#' Bootstrap uncertainty estimation for distance sampling models
#'
#' Performs a bootstrap for simple distance sampling models using the same data structures as \code{\link[mrds]{dht}}.
#'
#' @param model a model fitted by \code{\link{ds}} or a list of models
#' @param flatfile Data provided in the flatfile format. See \code{\link{flatfile}} for details.
#' @param convert.units conversion between units for abundance estimation, see "Units", below. (Defaults to 1, implying all of the units are "correct" already.) This takes precedence over any unit conversion stored in \code{model}.
#' @param resample_strata should resampling happen at the stratum (\code{Region.Label}) level? (Default \code{FALSE})
#' @param resample_obs should resampling happen at the observation (\code{object}) level? (Default \code{FALSE})
#' @param resample_transects should resampling happen at the transect (\code{Sample.Label}) level? (Default \code{TRUE})
#' @param nboot number of bootstrap replicates
#' @param summary_fun function that is used to obtain summary statistics from the bootstrap, see Summary Functions below. By default \code{\link{bootdht_Nhat_summarize}} is used, which just extracts abundance estimates.
#' @param select_adjustments select the number of adjustments in each bootstrap, when \code{FALSE} the exact detection function specified in \code{model} is fitted to each replicate. Setting this option to \code{TRUE} can significantly increase the runtime for the bootstrap. Note that for this to work \code{model} must have been fitted with \code{adjustment!=NULL}.
#' @param sample_fraction what proportion of the transects was covered (e.g., 0.5 for one-sided line transects).
#'
#' @section Summary Functions:
#' The function \code{summary_fun} allows the user to specify what summary statistics should be recorded from each bootstrap. The function should take two arguments, \code{ests} and \code{fit}. The former is the output from \code{dht2}, giving tables of estimates. The latter is the fitted detection function object. The function is called once fitting and estimation has been performed and should return a \code{data.frame}. Those \code{data.frame}s are then concatenated using \code{rbind}. One can make these functions return any information within those objects, for example abundance or density estimates or the AIC for each model. See Examples below.
#'
#' @section Model selection:
#' Model selection can be performed on a per-replicate basis within the bootstrap. This has three variations:
#' \enumerate{
#'    \item when \code{select_adjustments} is \code{TRUE} then adjustment terms are selected by AIC within each bootstrap replicate (provided that \code{model} had the \code{order} and \code{adjustment} options set to non-\code{NULL}.
#'    \item if \code{model} is a list of fitted detection functions, each of these is fitted to each replicate and results generated from the one with the lowest AIC.
#'    \item when \code{select_adjustments} is \code{TRUE} and \code{model} is a list of fitted detection functions, each model fitted to each replicate and number of adjustments is selected via AIC.
#' }
#' The last of these options can be very time consuming!
#' @importFrom utils txtProgressBar setTxtProgressBar getTxtProgressBar
#' @importFrom stats as.formula AIC
#' @importFrom mrds ddf dht
#' @seealso \code{\link{summary.dht_bootstrap}} for how to summarize the results, \code{\link{bootdht_Nhat_summarize}} for an example summary function.
#' @export
#' @examples
#' \dontrun{
#' # fit a model to the minke data
#' data(minke)
#' mod1 <- ds(minke)
#'
#' # summary function to save the abundance estimate
#' Nhat_summarize <- function(ests, fit) {
#'   return(data.frame(Nhat=ests$individuals$N$Estimate))
#' }
#'
#' # perform 5 bootstraps
#' bootout <- bootdht(mod1, flatfile=minke, summary_fun=Nhat_summarize, nboot=5)
#'
#' # obtain basic summary information
#' summary(bootout)
#' }
bootdht <- function(model,
                    flatfile,
                    resample_strata=FALSE,
                    resample_obs=FALSE,
                    resample_transects=TRUE,
                    nboot=100,
                    summary_fun=bootdht_Nhat_summarize,
                    convert.units=1,
                    select_adjustments=FALSE,
                    sample_fraction=1){

  if(!any(c(resample_strata, resample_obs, resample_transects))){
    stop("At least one of resample_strata, resample_obs, resample_transects must be TRUE")
  }

  # if we have a list...
  if(!all(class(model)=="list")){
    models <- list(model)
    # yes, I am a monster
  }else{
    models <- model
  }

  if(missing(convert.units)){
    convert.units <-  NULL
  }

  for(i in seq_along(models)){
    # only use valid ds models
    if(!all(class(models[[i]])=="dsmodel")){
      stop("Only models fitted by Distance::ds() may be used")
    }
  }
  dat <- flatfile

  # if we're using the default summary function and have Area 0 then
  # we're not going to have a good time
  if(missing(summary_fun) &
     (is.null(flatfile$Area) || all(flatfile$Area==0))){
    stop("No Area in flatfile, densities will be returned and the default summary function records only abundances. You need to write your own summary_fun.")
  }

  # apply the sample fraction
  check_sample_fraction(sample_fraction)
  dat$Effort <- dat$Effort*sample_fraction

  # this can be generalized later on
  stratum_label <- "Region.Label"
  obs_label <- "object"
  sample_label <- "Sample.Label"

  # which resamples are we going to make?
  possible_resamples <- c(stratum_label, sample_label, obs_label)
  our_resamples <- possible_resamples[c(resample_strata, resample_transects,
                                        resample_obs)]

  # count failures
  nbootfail <- 0
  # function to do a single bootstrap iteration
  bootit <- function(bootdat, our_resamples, groups,
                     convert.units, pb){

    # sample at the right levels
    for(sample_thingo in our_resamples){
      # what are the possible samples at this level
      levs <- unique(dat[[sample_thingo]])
      nlevs <- length(levs)
      levs <- sample(levs, nlevs, replace=TRUE)

      # make a new data frame with the correct number of replicates of the
      # per-stratum data in it
      bootdat <- lapply(levs, function(x) bootdat[bootdat[[sample_thingo]] == x, ])
      # make a special index to make unique IDs later
      iind <- rep(1:length(bootdat), lapply(bootdat, nrow))
      # make list of data.frames into one frame
      bootdat <- do.call("rbind", bootdat)
      # put that ID in there
      bootdat[[paste0(sample_thingo, "_ID")]] <- iind
    }

    # need unique object IDs
    bootdat$object <- 1:nrow(bootdat)
    # get the sample labels right
    bootdat$Sample.Label <- paste0(bootdat[[sample_label]], "-",
                                   bootdat[[paste0(sample_thingo, "_ID")]])


    aics <- rep(NA, length(models))
    for(i in seq_along(models)){
      model <- models[[i]]

      # setup the call to ds
      df_call <- model$call

      # if we want the number of adjustments to be selected each iteration...
      if(select_adjustments){
        df_call$order <- NULL
      }
      # insert the new data into the model
      df_call$data <- bootdat
      if(!is.null(convert.units)){
        df_call$convert.units <- convert.units
      }

      # fit that and update what's in models
      models[[i]] <- try(suppressMessages(eval(df_call, parent.frame(n=3))),
                         silent=TRUE)

      if(any(class(models[[i]]) == "try-error")){
        # if the model failed, return NA
        aics[i] <- NA
      }else{
        # if that wasn't bad, grab the AIC
        aics[i] <- AIC(models[[i]])$AIC
      }
    }

    # update progress bar
    setTxtProgressBar(pb, getTxtProgressBar(pb)+1)

    if(all(is.na(aics))){
      # if no models fitted, return NA
      nbootfail <<- nbootfail + 1
      return(NA)
    }else{
      fit <- models[[which.min(aics)]]
      # handle errors
      if(any(class(fit) == "try-error") ||
         any(is.na(fit$ddf$hessian))){
        nbootfail <<- nbootfail + 1
        return(NA)
      }else{
        return(summary_fun(fit$dht, fit$ddf))
      }
    }
  }

  cat(paste0("Performing ", nboot, " bootstraps\n"))
  pb <- txtProgressBar(0, nboot, style=3)
  # run the code
  boot_ests <- replicate(nboot,
                         bootit(dat, our_resamples,
                                summary_fun, convert.units=convert.units,
                                pb=pb), simplify=FALSE)

  # the above is then a list of thingos, do the "right" thing and assume
  # they are data.frames and then rbind them all together
  boot_ests <- do.call(rbind.data.frame, boot_ests)
  cat("\n")

  attr(boot_ests, "nboot") <- nboot
  attr(boot_ests, "nbootfail") <- nbootfail
  class(boot_ests) <- "dht_bootstrap"
  return(boot_ests)
}
