#' Convert units for abundance estimation
#'
#' It is often the case that effort, distances and prediction area are collected in different units in the field. Functions in \code{Distance} allow for an argument to convert between these and provide an answer that makes sense. This function calculates that conversion factor, given knowledge of the units of the quantities used.
#'
#' \code{convert_units} expects particular names for its inputs -- these should be singular names of the unit (e.g., "metre" rather than "metres"). You can view possible options with \code{\link{units_table}}. Both UK and US spellings are acceptable, case does not matter. For density estimation, area must still be provided ("objects per square ???"). Note that for cue counts (or other multiplier-based methods) one will still have to ensure that the rates are in the correct units for the survey.
#'
#'@param distance_units units distances were measured in.
#'@param effort_units units that effort were measured in. Set as \code{NULL} for point transects.
#'@param area_units units for the prediction area.
#'
#' @export
#' @author David L Miller
#'
#' @examples
#' # distances measured in metres, effort in kilometres and
#' # abundance over an area measured in hectares:
#' convert_units("Metre", "Kilometre", "Hectare")
#'
#' # all SI units, so the result is 1
#' convert_units("Metre", "metre", "square metre")
#'
#' # for points ignore effort
#' convert_units("Metre", NULL, "Hectare")
convert_units <- function(distance_units, effort_units, area_units){

  # get the unit table
  ut <- units_table()
  # make everything lower case to make matching easier
  ut$Unit <- tolower(ut$Unit)
  distance_units <- tolower(distance_units)
  area_units <- tolower(area_units)

  # grab the conversion factors from the table
  du <- ut$Conversion[ut$Unit==distance_units]
  au <- ut$Conversion[ut$Unit==area_units]

  # deal with point transects where effort is "number of visits"
  # but we have distance (truncation) squared
  if(is.null(effort_units)){
    eu <- ut$Conversion[ut$Unit==distance_units]
    pt <- TRUE
  }else{
    # otherwise do as above
    effort_units <- tolower(effort_units)
    eu <- ut$Conversion[ut$Unit==effort_units]
    pt <- FALSE
  }

  cu <- (du*eu)/au
  if(pt) cu <- sqrt(cu)

  # return the conversion
  return(cu)
}
