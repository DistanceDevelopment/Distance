#' Fit detection functions to line or point transect data
#'
#' BLURB
#'
#' @param data a \code{data.frame} containing at least a column called
#'        \code{distance}.
#' @param truncation truncation distance.
#' @param transect indicates transect type "line" (default) or "point".
#' @param formula formula for the scale parameter. For a CDS analysis leave 
#'        this as its default \code{~1}.
#' @param key key function to use; "hn" gives half-normal (default) or "hr"
#'        gives hazard-rate.
#' @param adjustment adjustment terms to use; "cos" gives cosine (default),
#'        "herm" gives Hermite polynomial and "poly" gives simple polynomial.
#'        "cos" is recommended.
#' @param order orders of the adjustment terms (as a vector), the default 
#'        value (\code{NULL}) will select via AIC.
#' TKTKTK talk about which orders are correct
#' @param scale the scale by which the distances in the adjustment terms are
#'        divided. Defaults to "scale", the scale parameter of the detection
#'        function. Other option is "width" which scales by the truncation
#'        distance.
#' @param binned if the data are binned, set this to \code{TRUE} (default 
#'        \code{FALSE}).
# @value
# TKTKTK need to think about this!!
#' @author David L. Miller
#'
#' @examples
#'  # An example from mrds, the golf tee data.
#'  library(distance)
#'  data(book.tee.data)
#'  tee.data<-book.tee.data$book.tee.dataframe
#'  ds.model<-ds(tee.data,4)
#'  summary(ds.model)
ds<-function(data, truncation, transect="line", formula=~1, key="hn",
             adjustment="cos", order=NULL, scale="scale",binned=FALSE){
  
  # this routine just creates a call to mrds, it's not very exciting
  # or fancy

  ## make sure that the data are in the right format first
  if(is.null(data$distance)){
    stop("Your data must have a column called 'distance'!\n")
  }

  # make sure that we have a data.frame()
  data<-data.frame(data)

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

  ## check inputs
  if(!(key %in% c("hn","hr"))){
    stop("key function must be \"hn\" or \"hr\".\n")
  }
  if(!(adjustment %in% c(NULL,"cos","herm","poly"))){
    stop("adjustment terms must be one of NULL, \"cos\", \"herm\" or \"poly\".\n")
  }
## adjustment order check here
## check for each adjustment type
## check set min.order for aic!
min.order <- 2
max.order <- 6
  if(is.null(truncation)){
    stop("Please supply truncation.\n")
  }

  # points or lines?
  if(transect=="line"){
    point <- FALSE
  }else if(transect=="point"){
    point <- TRUE
  }else{
    stop("Only \"point\" or \"line\" transects may be supplied.\n")
  }



  ### Actually fit some models here
  
  # if we are doing AIC search for adjustments
  aic.search <- FALSE
  if(!is.null(adjustment) & is.null(order)){
    aic.search <- TRUE
    order <- seq(min.order,max.order,by=2)
  }


  models<-list()
  aics<-c()

  for(i in seq(along=order)){ 

    ## CDS models
    if(as.character(formula)[2]=="1"){
      model.formula <- paste("~cds(key =\"", key,"\", formula = ~1,",
                            "adj.series=\"",adjustment,
                        "\",adj.order=",order[1:i],",",
                            "adj.scale=\"",scale,"\")",sep="")
    ## MCDS models
    }else{
      model.formula <- paste("~mcds(key = \"",key,"\",", 
                                   "formula =~",as.character(formula)[2],",", 
                                   "adj.series=\"",adjustment,"\",",
                                   "adj.order=",order[1:i],",",
                                   "adj.scale=\"",scale,"\")",sep="")
    }

    models[[i]]<-ddf(dsmodel = as.formula(model.formula),data = data, 
               method = "ds", meta.data = list(width = truncation,
                                               point = point,
                                               binned=binned))

    aics<-c(aics,models[[i]]$criterion)
  }

  model<-models[[which.min(aics)]]

  return(model)

}
