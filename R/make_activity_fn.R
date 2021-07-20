#' Multiplier bootstrap helper functions
#'
#' Use [activity::fitact] to fit an activity model.
#'
#' Uses [activity::fitact] to generate single possible availability estimates
#' based on bootstraps. The function returns another function, which can be
#' passed to `bootdht`. It is recommended that you try out the function before
#' passing it to `bootdht`. See examples for a template for use.
#'
#' @inheritParams activity::fitact
#' @param detector_daily_duration by default we assume that detectors were able to detect animals for 24 hours, if they were only able to do this for some proportion of the day (say daylight hours), then adjust this argument accordingly
#' @return a function which generates a single bootstrap estimate of
#' availability
#' @author David L Miller
#' @export
make_activity_fn <- function(..., detector_daily_duration=24){

  # check we can use activity
  if (!requireNamespace("activity", quietly = TRUE)){
    stop("Package 'activity' not installed!")
  }

  # save the arguments passed to the function
  call <- sys.call(sys.parent())
  # save detection duration separately, as it's not an arg to fitact
  detector_daily_duration <- call$detector_daily_duration
  call$detector_daily_duration <- NULL
  # actually do the call matching now
  args <- as.list(match.call(fitact, call=call, expand.dots=TRUE))
  # remove the function name
  args[[1]] <- NULL
  # ensure we only get one rep and don't show progress bar
  args$reps <- 1
  args$show <- FALSE

  # now eval all the args so that they are available to the parallel sessions
  for(i in seq(1, length(args))){
    args[[i]] <- eval(args[[i]], envir=parent.frame(n=2))
  }


  # function to return
  function(){
    # return the lower %ile for one rep, rescaling for proportion of day
    # where the detector was "on"
    unname(do.call(fitact, args)@act[3])/(detector_daily_duration/24)
  }
}
