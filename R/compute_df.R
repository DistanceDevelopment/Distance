#' @importFrom mrds DeltaMethod
compute_df <- function(k, type){
  df <- sapply(k, compute.df, type)
  df[df<1] <- 1
  return(df)
}

# Define function: compute.df
compute.df <- function(k, type){
  if(type=="O1" | type=="O2"| type=="O3"){
    H.O <- k - 1
    k.h.O <- rep(2, H.O)
    df <- sum(k.h.O - 1)
  }else{
    if(type=="S1" | type=="S2"){
      H.S <- floor(k/2)
      k.h.S <- rep(2, H.S)
      if(k %% 2 > 0) k.h.S[H.S] <- 3
      df <- sum(k.h.S - 1)
    }else{
      df <- k-1
    }
  }
  return(df)
}
