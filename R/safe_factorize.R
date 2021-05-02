# Safely make the right things factors in data.frames
#
# If you have a formula with as.factor(var) in it, this function will create a
# new column with the factor version of var with a column called
# as.factor(var). Uses model.frame
#
# @param formula an R formula
# @param data a data.frame
# @return data, with (maybe) some extra columns
# @author David L Miller
#' @importFrom stats model.frame
safe_factorize <- function(formula, data){

  if(!grepl("as\\.factor", as.character(formula)[2])){
    return(data)
  }

  if(class(try(model.frame(formula, data), silent=TRUE)) !=
       "try-error"){
    dn <- colnames(data)
    mf <- model.frame(formula, data)
    mf <- mf[, !(names(mf) %in% names(data)), drop=FALSE]
    data <- cbind.data.frame(data, mf)
  }

  return(data)
}
