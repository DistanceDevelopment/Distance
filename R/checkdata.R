#' Check that the data supplied to `ds` is correct
#'
#' This is an internal function that checks the `data.frame`s supplied
#' to `ds` are "correct".
#'
#' @param data as in `ds`
#' @param sample.table as in `ds`
#' @param region.table as in `ds`
#' @param obs.table as in `ds`
#' @param formula formula for the covariates
#'
#' @return Throws an error if something goes wrong, otherwise returns a
#' `data.frame`.
#'
#' @author David L. Miller
checkdata <- function(data, region.table=NULL, sample.table=NULL,
                      obs.table=NULL, formula=~1){

  # Check if the user has passed in a numeric vector
  if(!inherits(data, "data.frame") & !inherits(data, "list")){
    if(is.numeric(data)){
      data <- data.frame(distance = data)
    }else{
      stop("data is not a data.frame nor are the values numeric")
    }
  }

  # make sure that the data are in the right format first
  if(is.null(data$distance)){
    if(is.null(data$distend) & is.null(data$distbegin)){
      stop("Your data must (at least) have a column called 'distance' or 'distbegin' and 'distend'!")
    }else{
      data$distance <- (data$distend + data$distbegin)/2
    }
  }

  # make sure that we have a data.frame()
  data <- as.data.frame(data)

  # first see if the data has detected/observer/object columns, if not add them
  if(!any(names(data) == "observer")){
    data <- cbind(data, observer=rep(1, nrow(data)))
  }else if(!all(data$observer %in% c(1, 2))){
    stop("'observer' is a reserved column name, please rename this column.")
  }
  if(!any(names(data) == "detected")){
    data <- cbind(data, detected=rep(1, nrow(data)))
  }
  if(!any(names(data) == "object")){
    data <- cbind(data, object=1:nrow(data))
  }else{
    # check that the object IDs are unique
    # first need to remove the rows with NA distances used for padding
    # below
    data_no_NA <- data[!is.na(data$distance), ]
    if(length(data_no_NA$object) != length(unique(data_no_NA$object))){
      stop("Not all object IDs are unique, check data.")
    }
  }

  ## check that the covariates that are specified exist in the data
  formula <- as.formula(formula)
  if(formula!=~1){
    vars_in_data <- all.vars(formula) %in% names(data)
    if(!all(vars_in_data)){
      stop(paste("Variable(s):",
                 paste(all.vars(formula)[!vars_in_data],collapse=", "),
                 "are in the model formula but not in the data."))
    }
  }

  if(is.null(region.table) & is.null(sample.table) & is.null(obs.table)){
    ## if the tables are NULL then we either have a detection function only
    ## or we have a simplified table, in which case we need to interpret into
    ## three tables.

    if(all(c("Region.Label", "Sample.Label", "Effort", "object") %in%
           colnames(data))){
      ## construct region table
      # if area isn't specified, we can still compute density
      if(!("Area" %in% names(data))){
        message("No Area column found, generating density estimates only")
        data$Area <- 0
      }
      region.table <- unique(data[,c("Region.Label", "Area")])
      # make sure that the region areas are consistent -- the above can
      # lead to duplicate labels if the areas are not the same for a given
      # region label
      if(nrow(region.table) != length(unique(data$Region.Label))){
        stop("Region areas are not consistent.")
      }
      rownames(region.table) <- 1:nrow(region.table)
      # drop Area column
      data <- data[, !c(colnames(data) %in% "Area")]

      ## construct sample table
      sample.table <- unique(data[, c("Sample.Label", "Region.Label",
                                      "Effort")])
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
      data <- data[, !c(colnames(data) %in% "Effort")]


      ## construct obs
      obs.table <- unique(data[, c("object", "Region.Label","Sample.Label")])
      # remove obs/region combinations with no samples in them
      obs.table <- obs.table[!is.na(obs.table$Sample.Label), ]
      rownames(obs.table) <- 1:nrow(obs.table)

      # drop Region and Sample label columns
      # actually don't do this because then we can't use subset= in dht
      #data <- data[,!c(colnames(data) %in% c("Region.Label","Sample.Label"))]
      rownames(data) <- 1:nrow(data)

      # remove the NA rows
      data <- data[!is.na(data$distance),]
    }else if(all(tolower(c("Region.Label", "Sample.Label",
                           "Effort", "object")) %in%
                 tolower(colnames(data)))){
      stop("flatfile column names detected but with incorrect cAsE, correct the names and re-run. See ?flatfile for details.")
    }
  }else{
    # check that dht info has the right column titles
    if(!is.null(region.table) & !is.null(sample.table) & !is.null(obs.table)){
      if(!("Area" %in% names(region.table))){
        message("No Area column found, generating density estimates only")
        region.table$Area <- 0
      }
      if(!all(c("Region.Label", "Area") %in% names(region.table))){
        stop("region.table must have columns named 'Region.Label' and 'Area'")
      }
      if(!all(c("Region.Label", "Sample.Label", "Effort") %in%
              names(sample.table))){
        stop("sample.table must have columns named 'Region.Label', 'Sample.Label' and 'Effort'")
      }
      if(!all(c("Region.Label", "Sample.Label", "object") %in%
              names(obs.table))){
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
