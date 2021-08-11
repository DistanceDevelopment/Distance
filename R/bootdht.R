#' Bootstrap uncertainty estimation for distance sampling models
#'
#' Performs a bootstrap for simple distance sampling models using the same data
#' structures as [`dht`][mrds::dht]. Note that only geographical stratification
#' as supported in `dht` is allowed.
#'
#' @param model a model fitted by [`ds`] or a list of models
#' @param flatfile Data provided in the flatfile format. See [`flatfile`] for
#' details.
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
#' [`bootdht_Nhat_summarize`] is used, which just extracts abundance estimates.
#' @param select_adjustments select the number of adjustments in each
#' bootstrap, when `FALSE` the exact detection function specified in `model` is
#' fitted to each replicate. Setting this option to `TRUE` can significantly
#' increase the runtime for the bootstrap. Note that for this to work `model`
#' must have been fitted with `adjustment!=NULL`.
#' @param sample_fraction what proportion of the transects was covered (e.g.,
#' 0.5 for one-sided line transects).
#' @param multipliers `list` of multipliers. See "Multipliers" below.
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
#' @section Multipliers:
#' It is often the case that we cannot measure distances to individuals or
#' groups directly, but instead need to estimate distances to something they
#' produce (e.g., for whales, their blows; for elephants their dung) -- this is
#' referred to as indirect sampling. We may need to use estimates of production
#' rate and decay rate for these estimates (in the case of dung or nests) or
#' just production rates (in the case of songbird calls or whale blows). We
#' refer to these conversions between "number of cues" and "number of animals"
#' as "multipliers".
#'
#' The `multipliers` argument is a `list`, with 3 possible elements (`creation`
#' and `decay`). Each element of which is either:
#' * `data.frame` and must have at least a column named `rate`, which abundance
#'    estimates will be divided by (the term "multiplier" is a misnomer, but
#'    kept for compatibility with Distance for Windows). Additional columns can
#'    be added to give the standard error and degrees of freedom for the rate
#'    if known as `SE` and `df`, respectively. You can use a multirow
#'    `data.frame` to have different rates for different geographical areas
#'    (for example). In this case the rows need to have a column (or columns)
#'    to `merge` with the data (for example `Region.Label`).
#' * a `function` which will return a single estimate of the relevant
#'   multiplier. See [`make_activity_fn`] for a helper function for use with the
#'   `activity` package.
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
#' If `cores`>1 then the `parallel`/`doParallel`/`foreach`/`doRNG` packages
#' will be used to run the computation over multiple cores of the computer. To
#' use this component you need to install those packages using:
#' `install.packages(c("foreach", "doParallel", "doRNG"))` It is advised that
#' you do not set `cores` to be greater than one less than the number of cores
#' on your machine. The `doRNG` package is required to make analyses
#' reproducible ([`set.seed`] can be used to ensure the same answers).
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
#' @seealso [`summary.dht_bootstrap`] for how to summarize the results,
#' [`bootdht_Nhat_summarize`] for an example summary function.
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
                    multipliers=NULL,
                    progress_bar="base",
                    cores=1){

  if(!any(c(resample_strata, resample_obs, resample_transects))){
    stop("At least one of resample_strata, resample_obs, resample_transects must be TRUE")
  }

  # make everything a list...
  if(!all(class(model)=="list")){
    models <- list(model)
    # yes, I am a monster
  }else{
    models <- model
  }

  if(missing(convert.units)){
    convert.units <-  NULL
  }

  # only use valid ds models
  for(i in seq_along(models)){
    if(!all(class(models[[i]])=="dsmodel")){
      stop("Only models fitted by Distance::ds() may be used")
    }
  }
  dat <- flatfile
  if(!("object" %in% names(dat))){
    dat$object <- 1:nrow(dat)
  }

  # if we're using the default summary function and have Area 0 then
  # we're not going to have a good time
  if(missing(summary_fun) &
     (is.null(flatfile$Area) || all(flatfile$Area==0))){
    stop("No Area in flatfile, densities will be returned and the default summary function records only abundances. You need to write your own summary_fun.")
  }

  # apply the sample fraction
  dat <- dht2_sample_fraction(sample_fraction, dat)
  dat$Effort <- dat$Effort * dat$sample_fraction
  dat$sample_fraction <- NULL

  # process non-function multipliers
  multipliers_nofun <- Filter(Negate(is.function), multipliers)
  if(length(multipliers_nofun) > 0){
    dat <- dht2_multipliers(multipliers_nofun, dat)
  }else{
    dat$rate <- 1
  }

  # get any multiplier functions
  multipliers_fun <- Filter(is.function, multipliers)

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
        model$call[[names(lpars)[i]]] <- eval(lpars[[i]],
                                              envir=parent.frame(n=3))
      }
    }
    model
  })

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

  # run the bootstrap
  if(cores > 1){
    if (!requireNamespace("foreach", quietly = TRUE) &
        !requireNamespace("doParallel", quietly = TRUE) &
        !requireNamespace("doRNG", quietly = TRUE) &
        !requireNamespace("parallel", quietly = TRUE)){
      stop("Packages 'parallel', 'foreach' and 'doParallel' need to be installed to use multiple cores.")
    }

    # build the cluster
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    # shutdown cluster when the function exits
    # (this works even if the function crashes)
    on.exit(parallel::stopCluster(cl))

    # load the activity package in the other sessions
    if(length(multipliers_fun) > 0){
      packages <- c("activity")
    }else{
      packages <- NULL
    }

    # needed to avoid a syntax error/check fail
    `%dopar%` <- foreach::`%dopar%`
    `%dorng2%` <- doRNG::`%dorng%`
    # fit the model nboot times over cores cores
    # note there is a bit of fiddling here with the progress bar to get it to
    # work (updates happen in this loop rather than in bootit)
    boot_ests <- foreach::foreach(i=1:nboot, .packages=packages) %dorng2% {
      bootit(dat, models=models, our_resamples=our_resamples,
             summary_fun=summary_fun, convert.units=convert.units,
             pb=list(increment=function(pb){invisible()}),
             multipliers_fun=multipliers_fun, sample_label=sample_label,
             select_adjustments=select_adjustments)
    }

  }else{
    boot_ests <- replicate(nboot,
                           bootit(dat, models=models, our_resamples,
                                  summary_fun, convert.units=convert.units,
                                  pb=pb, multipliers_fun=multipliers_fun,
                                  sample_label=sample_label,
                                  select_adjustments=select_adjustments),
                           simplify=FALSE)
  }

  # do some post-processing
  fail_fun <- function(x) class(x)=="bootstrap_failure"
  # add replicate IDs
  bootids <- seq_len(length(boot_ests))
  # how many failures
  failures <- length(Filter(fail_fun, boot_ests))
  # remove failures from IDs
  bootids <- as.list(bootids[unlist(lapply(boot_ests, Negate(fail_fun)))])
  # get just the "good" results
  boot_ests <- Filter(Negate(fail_fun), boot_ests)
  # add IDs to estimates list
  boot_ests <- mapply(cbind.data.frame,
                      boot_ests, bootstrap_ID=bootids,
                      SIMPLIFY=FALSE)

  # the above is then a list of thingos, do the "right" thing and assume
  # they are data.frames and then rbind them all together
  boot_ests <- do.call(rbind.data.frame, boot_ests)

  cat("\n")

  attr(boot_ests, "nboot") <- nboot
  attr(boot_ests, "failures") <- failures
  class(boot_ests) <- c("dht_bootstrap", "data.frame")
  return(boot_ests)
}
