get_truncation <- function(truncation, cutpoints, data){
  
  # Check truncation is specified
  if (is.null(truncation)) {
    stop("Please supply truncation distance or percentage.", call. = FALSE)
  }
  
  # Check the structure of the truncation
  if (!(length(truncation) == 1 || (is.list(truncation) &&
                                    length(truncation) == 2 &&
                                    all(c("left", "right") %in% names(truncation)) &&
                                    all(sapply(truncation, function(el) length(el) == 1))))) {
    stop("Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\", each being a single value.", call. = FALSE)
  }
  
  # Check if truncation is character (i.e., contains %) and if used with binned data
  if (any(sapply(truncation, is.character)) &&
      (!is.null(cutpoints) || any(c("distbegin", "distend") %in% colnames(data)))) {
    stop("Truncation cannot be supplied as a percentage with binned data.", call. = FALSE)
  }
  
  # Check if any truncation value is a character but doesn't contains the "%" symbol
  if(any(sapply(truncation, function(x){is.character(x) && !grepl("%", x)}))){
    warning("Truncation values supplied as characters will be interpreted as % truncation values.", call. = FALSE)
  }
  
  # Helper function to convert percentage strings to quantiles
  percent_trunc <- function(value, side = c("left", "right")) {
    value_num <- as.numeric(sub("%", "", value))
    prob <- if (side == "left") value_num / 100 else 1 - (value_num / 100)
    quantile(data$distance, probs = prob, na.rm = TRUE)
  }
  
  # Get the truncation values
  if(is.list(truncation)){
    
    # Left truncation
    if(is.numeric(truncation$left)){
      left <- truncation$left
    } else if(is.character(truncation$left)){
      left <- percent_trunc(truncation$left, "left")
    } else stop("Left truncation must be a number or percentage string.")
    
    # Right truncation
    if(is.numeric(truncation$right)){
      width <- truncation$right
    } else if(is.character(truncation$right)){
      width <- percent_trunc(truncation$right, "right")
    } else stop("Right truncation must be a number or percentage string.")
    
  #If it is not a list
  } else if(is.numeric(truncation)) {
    width <- truncation
    left <- NULL
  } else if (is.character(truncation)){
    width <- percent_trunc(truncation, "right")
    left <- NULL
  } else stop("Right truncation must be a number or percentage string.")
  
  # Now check that truncation is not greater than largest cutpoint if binned data
  if(!is.null(cutpoints)){
    if(width > cutpoints[length(cutpoints)]){
      warning(paste("Truncation width is greater than the largest bin distance, re-setting truncation to be largest cutpoint value: ", cutpoints[length(cutpoints)], sep = ""), immediate. = TRUE, call. = FALSE)
      # Make truncation largest cutpoint
      width <- cutpoints[length(cutpoints)]
    }
  }
  
  list(left=left, width=width)
}