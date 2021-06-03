get_truncation <- function(truncation, cutpoints, data){

  if(is.null(truncation)){
    stop("Please supply truncation distance or percentage.")
  }else if(any(unlist(lapply(truncation, is.character))) &
           (!is.null(cutpoints) |
            any(c("distbegin", "distend") %in% colnames(data))
           )){
    stop("Truncation cannot be supplied as a percentage with binned data")
  }else{
    # if we have left truncation too...
    if(is.list(truncation)){
      if((any(names(truncation)=="left") &
          any(names(truncation)=="right")) &
         length(truncation)==2){

        # check for each of left and right that we have % or distance...
        # left
        if(is.double(truncation$left) & length(truncation$left)==1){
          left <- truncation$left
        }else if(is.character(truncation$left) & length(truncation$left)==1){
          # % string to number
          truncation$left <- as.numeric(sub("%","",truncation$left))
          left <- quantile(data$distance, probs=truncation$left/100,
                           na.rm=TRUE)
        }else{
          stop("Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
        }
        # right
        if(is.double(truncation$right) & length(truncation$right)==1){
          width <- truncation$right
        }else if(is.character(truncation$right) & length(truncation$right)==1){
          # % string to number
          truncation$right <- as.numeric(sub("%", "", truncation$right))
          width <- quantile(data$distance, probs=1-(truncation$right/100),
                            na.rm=TRUE)
        }else{
          stop("Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
        }
      }else{
        stop("Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
      }

    # just right truncation
    }else if(is.numeric(truncation) & length(truncation)==1){
      width <- truncation
      left <- NULL
    }else if(is.character(truncation) & length(truncation)==1){
      # % string to number
      truncation <- as.numeric(sub("%","",truncation))
      width <- quantile(data$distance, probs=1-(truncation/100), na.rm=TRUE)
      left <- NULL
    }else{
      stop("Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
    }
  }
  list(left=left, width=width)
}
