# temporary internal function to deal with deprecated arguments
.deprecated_args <- function(deprecated, args){

  bad <- names(args) %in% deprecated

  if(sum(bad)==1){
    stop(paste0("Argument: ",
                paste(names(args)[bad], collapse=", "),
                " is deprecated, check documentation."))
  }else if(sum(bad)>1){
    stop(paste0("Argument(s): ",
                paste(names(args)[bad], collapse=", "),
                " are deprecated, check documentation."))
  }

}
