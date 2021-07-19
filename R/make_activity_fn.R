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
#' @return a function which generates a single bootstrap estimate of
#' availability
#' @author David L Miller
#' @export
make_activity_fn <- function(...){

  # check we can use activity
  if (!requireNamespace("activity", quietly = TRUE)){
    stop("Package 'activity' not installed!")
  }

  # save the arguments passed to the function
  args <- as.list(match.call(fitact, expand.dots=TRUE))
  # remove the function name
  args[[1]] <- NULL
  # ensure we only get one rep and don't show progress bar
  args$reps <- 1
  args$show <- FALSE

  # function to return
  function(){
    # return the lower %ile for one rep
    unname(do.call(fitact, args)@act[3])
  }
}
