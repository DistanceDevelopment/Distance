#' Bootstrap uncertainty estimation for distance sampling models
#'
#' Performs a bootstrap for simple distance sampling models using the same data
#' structures as [`dht`][mrds::dht]. Note that only geographical stratification
#' as supported in `dht` is allowed.
#'
#' @param model a model fitted by [`ds`][ds] or a list of models
#' @param flatfile Data provided in the flatfile format. See
#' [`flatfile`][flatfile] for details.
#' @param convert.units conversion between units for abundance estimation, see
#' "Units", below. (Defaults to 1, implying all of the units are "correct"
#' already.) This takes precedence over any unit conversion stored in `model`.
#' @param resample_strata should resampling happen at the stratum
#' (`Region.Label`) level? (Default `FALSE`)
#' @param resample_obs should resampling happen at the observation (`object`)
#' level? (Default `FALSE`)
#' @param resample_transects should resampling happen at the transect
#' (`Sample.Label`) level? (Default `TRUE`)
#' @param nboot number of bootstrap replicates
#' @param summary_fun function that is used to obtain summary statistics from
#' the bootstrap, see Summary Functions below. By default
#' [`bootdht_Nhat_summarize`][bootdht_Nhat_summarize] is used, which just
#' extracts abundance estimates.
#' @param select_adjustments select the number of adjustments in each
#' bootstrap, when `FALSE` the exact detection function specified in `model` is
#' fitted to each replicate. Setting this option to `TRUE` can significantly
#' increase the runtime for the bootstrap. Note that for this to work `model`
#' must have been fitted with `adjustment!=NULL`.
#' @param sample_fraction what proportion of the transects was covered (e.g.,
#' 0.5 for one-sided line transects).
#' @param progress_bar which progress bar should be used? Default "base" uses
#' `txtProgressBar`, "none" suppresses output, "progress" uses the
#' `progress` package, if installed.
#' @param cores number of CPU cores to use to compute the estimates. See "Parallelization" below.
#'
#' @section Summary Functions:
#' The function `summary_fun` allows the user to specify what summary
#' statistics should be recorded from each bootstrap. The function should take
#' two arguments, `ests` and `fit`. The former is the output from
#' `dht2`, giving tables of estimates. The latter is the fitted detection
#' function object. The function is called once fitting and estimation has been
#' performed and should return a `data.frame`. Those `data.frame`s
#' are then concatenated using `rbind`. One can make these functions
#' return any information within those objects, for example abundance or
#' density estimates or the AIC for each model. See Examples below.
#'
#' @section Model selection:
#' Model selection can be performed on a per-replicate basis within the
#' bootstrap. This has three variations:
#'   1. when `select_adjustments` is `TRUE` then adjustment terms are selected
#'   by AIC within each bootstrap replicate (provided that `model` had the
#'   `order` and `adjustment` options set to non-`NULL`.
#'   2. if `model` is a list of fitted detection functions, each of these is
#'   fitted to each replicate and results generated from the one with the
#'   lowest AIC.
#'   3. when `select_adjustments` is `TRUE` and `model` is a list of fitted
#'   detection functions, each model fitted to each replicate and number of
#'   adjustments is selected via AIC.
#' This last option can be extremely time consuming.
#'
#' @section Parallelization:
#' If `cores`>1 then the `parallel`/`doParallel`/`foreach` packages will be
#' used to run the computation over multiple cores of the computer. To use this
#' component you need to install those packages using:
#' `install.packages(c("foreach", "doParallel"))` It is advised that you do not
#' set `cores` to be greater than one less than the number of cores on your
#' machine.
#'
#' It is also hard to debug any issues in `summary_fun` so it is best to run a
#' small number of bootstraps first in parallel to check that things work. On
#' Windows systems `summary_fun` does not have access to the global environment
#' when running in parallel, so all computations must be made using only its
#' `ests` and `fit` arguments (i.e., you can not use R objects from elsewhere
#' in that function, even if they are available to you from the console).
#'
#' @importFrom utils txtProgressBar setTxtProgressBar getTxtProgressBar
#' @importFrom stats as.formula AIC
#' @importFrom mrds ddf dht
#' @seealso [`summary.dht_bootstrap`][summary.dht_bootstrap] for how to
#' summarize the results, [`bootdht_Nhat_summarize`][bootdht_Nhat_summarize]
#' for an example summary function.
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
                    sample_fraction=1,
                    progress_bar="base",
                    cores=1){

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
  dat <- dht2_sample_fraction(sample_fraction, dat)

  # this can be generalized later on
  stratum_label <- "Region.Label"
  obs_label <- "object"
  sample_label <- "Sample.Label"

  # which resamples are we going to make?
  possible_resamples <- c(stratum_label, sample_label, obs_label)
  our_resamples <- possible_resamples[c(resample_strata, resample_transects,
                                        resample_obs)]

  # process models
  # this resolves all symbols in the call so arguments can be accessed when
  # running in parallel
  models <- lapply(models, function(model){
    lpars <- as.list(model$call)
    for(i in seq(2, length(lpars))){
      if(is.symbol(model$call[[names(lpars)[i]]])){
        model$call[[names(lpars)[i]]] <- eval(lpars[[i]])
      }
    }
    model
  })

  # count failures
  nbootfail <- 0
  # function to do a single bootstrap iteration
  bootit <- function(bootdat, models, our_resamples, summary_fun,
                     convert.units, pb, ...){
    # sample at the right levels
    for(sample_thingo in our_resamples){
      # what are the possible samples at this level
      levs <- unique(dat[[sample_thingo]])
      nlevs <- length(levs)
      levs <- sample(levs, nlevs, replace=TRUE)

      # make a new data frame with the correct number of replicates of the
      # per-stratum data in it
      bootdat <- lapply(levs, function(x){
        bootdat[bootdat[[sample_thingo]] == x, ]
      })
      # make a special index to make unique IDs later
      iind <- rep(seq_len(length(bootdat)), lapply(bootdat, nrow))
      # make list of data.frames into one frame
      bootdat <- do.call("rbind", bootdat)
      # put that ID in there
      bootdat[[paste0(sample_thingo, "_ID")]] <- iind
    }

    # need unique object IDs
    bootdat$object <- seq_len(nrow(bootdat))
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
      models[[i]] <- try(suppressMessages(eval(df_call)),
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
    pb$increment(pb$pb)

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

  if(cores > 1 & progress_bar != "none"){
    progress_bar <- "none"
    message("Progress bars cannot be shown when using cores>1")
  }

  # decide how to report progress
  if(progress_bar == "base"){
    pb <- list(pb        = txtProgressBar(0, nboot, style=3),
               increment = function(pb){
                 setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
               },
               set = function(pb, n, max) setTxtProgressBar(pb, n),
               done      = function(pb){
                 setTxtProgressBar(pb, environment(pb$up)$max)
               })
  }else if(progress_bar == "none"){
    pb <- list(pb        = NA,
               set = function(pb, n, max) invisible(),
               increment = function(pb) invisible(),
               done = function(pb) invisible())
  }else if(progress_bar == "progress"){
    if (!requireNamespace("progress", quietly = TRUE)){
      stop("Package 'progress' not installed!")
    }else{
      pb <- list(pb = progress::progress_bar$new(
                                       format="   [:bar] :percent eta: :eta",
                                       total=nboot, clear=FALSE, width=80),
                 set = function(pb, n, max) pb$update(n/max),
                 increment = function(pb) pb$tick(),
                 done = function(pb) pb$update(1))
      pb$pb$tick(0)
    }
  }else{
    stop("progress_bar must be one of \"none\", \"base\" or \"progress\"")
  }

  # run the code
  if(cores > 1){
    if (!requireNamespace("foreach", quietly = TRUE) &
        !requireNamespace("doParallel", quietly = TRUE) &
        !requireNamespace("parallel", quietly = TRUE)){
      stop("Packages 'parallel', 'foreach' and 'doParallel' need to be installed to use multiple cores.")
    }

    # build the cluster
    cl <- parallel::makeCluster(cores, outfile="")
    doParallel::registerDoParallel(cl)

    # needed to avoid a syntax error/check fail
    `%dopar2%` <- foreach::`%dopar%`
    # fit the model nboot times over cores cores
    # note there is a bit of fiddling here with the progress bar to get it to
    # work (updates happen in this loop rather than in bootit)
    boot_ests <- foreach::foreach(i=1:nboot) %dopar2% {
      bootit(dat, models=models, our_resamples=our_resamples,
             summary_fun=summary_fun, convert.units=convert.units,
             pb=list(increment=function(pb){invisible()}))
    }
    # shutdown cluster
    parallel::stopCluster(cl)

    # post-process
    nbootfail <- sum(unlist(lapply(boot_ests, is.na)))
    boot_ests <- Filter(function(x) !is.na(x), boot_ests)

  }else{
    boot_ests <- replicate(nboot,
                           bootit(dat, models=models, our_resamples,
                                  summary_fun, convert.units=convert.units,
                                  pb=pb), simplify=FALSE)
  }
  # the above is then a list of thingos, do the "right" thing and assume
  # they are data.frames and then rbind them all together
  boot_ests <- do.call(rbind.data.frame, boot_ests)

  cat("\n")

  attr(boot_ests, "nboot") <- nboot
  attr(boot_ests, "nbootfail") <- nbootfail
  class(boot_ests) <- "dht_bootstrap"
  return(boot_ests)
}
