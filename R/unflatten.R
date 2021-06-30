#' Unflatten flatfile data.frames
#'
#' Sometimes data is provided in the [`flatfile`][flatfile] format, but we
#' really want it in `mrds` format (that is, as distance data, observation
#' table, sample table and region table format). This function undoes the
#' flattening, assuming that the data have the correct columns.
#'
#' @param data data in flatfile format (a `data.frame`)
#' @return `list` of four `data.frame`s: distance data, observation table,
#' sample table, region table.
#'
#' @author David L Miller
#' @export
unflatten <- function(data){

  ## construct region table
  region.table <- unique(data[,c("Region.Label", "Area")])
  # make sure that the region areas are consistent -- the above can
  # lead to duplicate labels if the areas are not the same for a given
  # region label
  if(nrow(region.table) != length(unique(data$Region.Label))){
    stop("Region areas are not consistent.")
  }
  rownames(region.table) <- 1:nrow(region.table)
  # drop Area column
  data <- data[,!c(colnames(data) %in% "Area")]

  ## construct sample table
  sample.table <- unique(data[,c("Sample.Label", "Region.Label", "Effort")])
  # remove sample/region combinations with no samples in them
  sample.table <- sample.table[!is.na(sample.table$Effort), ]

  # possible that Effort is not the same for a given
  # Sample.Label+Region.Label -- this is BAD.
  if(nrow(sample.table)!=nrow(unique(sample.table[,c("Sample.Label",
                                              "Region.Label")]))){
    stop("A sample has a non-unique \"Effort\", check data!")
  }

  rownames(sample.table) <- 1:nrow(sample.table)
  # drop Effort column
  data <- data[,!c(colnames(data) %in% "Effort")]


  ## construct obs
  obs.table <- unique(data[,c("object", "Region.Label", "Sample.Label")])
  rownames(obs.table) <- 1:nrow(obs.table)

  # drop Region and Sample label columns
  # actually don't do this because then we can't use subset= in dht
  #data <- data[,!c(colnames(data) %in% c("Region.Label","Sample.Label"))]
  rownames(data) <- 1:nrow(data)

  # remove the NA rows
  data <- data[!is.na(data$distance), ]

  return(list(data=data, sample.table=sample.table, region.table=region.table,
              obs.table=obs.table))
}
