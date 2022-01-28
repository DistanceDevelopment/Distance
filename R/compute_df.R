# note this differs from the compute.df function in mrds
# as this version is vectorised over k and type, which should
# be the same length
compute_df <- function(k, type){

  if(length(k) != length(type)){
    stop("k and type must be the same length")
  }

  # storage
  df <- rep(NA, length(k))

  # handle O* estimators
  o_ind <- type=="O1" | type=="O2"| type=="O3"
  if(any(o_ind)){
    H.O <- k[o_ind] - 1
    k.h.O <- rep(2, H.O)
    df[o_ind] <- sum(k.h.O - 1)
  }

  # handle S* estimators
  s_ind <- type=="S1" | type=="S2"
  if(any(s_ind)){
    H.S <- floor(k[s_ind]/2)
    k.h.S <- rep(2, H.S)
    if(k[s_ind] %% 2 > 0) k.h.S[H.S] <- 3
    df[s_ind] <- sum(k.h.S - 1)
  }

  # everything else
  df[is.na(df)] <- k[is.na(df)]-1

  # must have at least one degree of freedom!
  df[df<1] <- 1
  return(df)
}
