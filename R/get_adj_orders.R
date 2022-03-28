# generate the appropriate orders for adjustments
get_adj_orders <- function(n, key, adjustment){

  # this is according to p. 47 of IDS
  if(key=="unif" & adjustment=="cos"){
    # for Fourier...
    order <- 1:n
  }else if(adjustment %in% c("poly", "herm")){
    # simple and Hermite poly: even from 4 (unless uniform)
    starto <- 4
    if(key == "unif") starto <- 2
    order <- seq(starto, by=2, length.out=n)
  }else if(adjustment=="cos"){
    # cosine: by 1 from 2
    order <- seq(2, by=1, length.out=n)
  }else{
    stop("Bad adjustment term definition")
  }

  return(order)
}
