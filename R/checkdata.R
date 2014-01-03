#' Check that the data supplied to \code{ds} is correct
#'
#' This is an internal function that checks the \code{data.frame}s supplied
#' to \code{ds} are "correct".
#'
#' @param data as in \code{ds}
#' @param sample.table as in \code{ds}
#' @param region.table as in \code{ds}
#' @param obs.table as in \code{ds}
#'
#' @return Throws an error if something goes wrong, otherwise returns a
#'  \code{data.frame}.
#'
#' @author David L. Miller
checkdata<-function(data,region.table=NULL,sample.table=NULL,obs.table=NULL){

   ## make sure that the data are in the right format first
  if(is.null(data$distance)){
    stop("Your data must (at least) have a column called 'distance'!")
  }

  # make sure that we have a data.frame()
  data <- as.data.frame(data)

  # first see if the data has detected/observer/object columns, if not add them
  if(!any(names(data)=="observer")){
    data<-cbind(data,observer=rep(1,nrow(data)))
  }
  if(!any(names(data)=="detected")){
    data<-cbind(data,detected=rep(1,nrow(data)))
  }
  if(!any(names(data)=="object")){
    data<-cbind(data,object=1:nrow(data))
  }

  if(is.null(region.table) & is.null(sample.table) & is.null(obs.table)){
    ## if the tables are NULL then we either have a detection function only
    ## or we have a simplified table, in which case we need to interpret into
    ## three tables.

    if(all(c("Region.Label", "Area", "Sample.Label", "Effort", "object") %in%
           colnames(data))){
      ## construct region table
      region.table <- unique(data[,c("Region.Label", "Area")])
      rownames(region.table) <- 1:nrow(region.table)
      # drop Area column
      data <- data[,!c(colnames(data) %in% "Area")]

      ## construct sample table
      sample.table <- unique(data[,c("Sample.Label", "Region.Label", "Effort")])
      rownames(sample.table) <- 1:nrow(sample.table)
      # drop Effort column
      data <- data[,!c(colnames(data) %in% "Effort")]


      ## construct obs
      obs.table <- unique(data[,c("object", "Region.Label","Sample.Label")])
      rownames(obs.table) <- 1:nrow(obs.table)
      # drop Region and Sample label columns
      data <- data[,!c(colnames(data) %in% c("Region.Label","Sample.Label"))]
      rownames(data) <- 1:nrow(data)
    }

  }else{
    # check that dht info has the right column titles
    if(!is.null(region.table) & !is.null(sample.table) & !is.null(obs.table)){
      if(!all(c("Region.Label","Area") %in% names(region.table))){
        stop("region.table must have columns named 'Region.Label' and 'Area'")
      }
      if(!all(c("Region.Label","Sample.Label","Effort") %in%
              names(sample.table))){
        stop("sample.table must have columns named 'Region.Label', 'Sample.Label' and 'Effort'")
      }
      if(!all(c("Region.Label","Sample.Label","object") %in% names(obs.table))){
        stop("obs.table must have columns names 'Region.Label', 'Sample.Label' and 'object'")
      }
    }
  }

  # nothing bad happened, yay!
  return(list(data         = data,
              region.table = region.table,
              sample.table = sample.table,
              obs.table    = obs.table))
}
