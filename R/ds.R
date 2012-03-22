#' Fit detection functions to line or point transect data
#'
#' Some blurb here about distance sampling analysis. 
#'
#' @param data a \code{data.frame} containing at least a column called
#'        \code{distance}.
#' @param truncation either truncation distance (numeric, e.g. 5) or percentage
#'        (as a string, e.g. "15\%"). Can be supplied as a \code{list}
#'        with elements \code{left} and \code{right} if left truncation is
#'        required (e.g. \code{list(left=1,right=20)} or 
#'        \code{list(left="1\%",right="15\%")} or even 
#         \code{list(left="1",right="15\%")}).
#'        When specified as a percentage, the largest \code{right} and 
#'        smallest \code{left} percent distances are discarded. 
#' @param transect indicates transect type "line" (default) or "point".
#' @param formula formula for the scale parameter. For a CDS analysis leave 
#'        this as its default \code{~1}.
#' @param key key function to use; "hn" gives half-normal (default) or "hr"
#'        gives hazard-rate.
#' @param adjustment adjustment terms to use; "cos" gives cosine (default),
#'        "herm" gives Hermite polynomial and "poly" gives simple polynomial.
#'        "cos" is recommended.
#' @param order orders of the adjustment terms to fit (as a vector/scalar), the
#'        default value (\code{NULL}) will select via AIC. For cosine 
#'        adjustments, valid orders are integers greater than 2. For Hermite 
#'        polynomials, even integers equal or greater than 4 are allowed. For 
#'        simple polynomials even integers equal or greater than 2 are allowed.
#' @param scale the scale by which the distances in the adjustment terms are
#'        divided. Defaults to "scale", the scale parameter of the detection
#'        function. Other option is "width" which scales by the truncation
#'        distance.
#' @param cutpoints if the data are binned, this vector gives the cutpoints of 
#'        the bins. Ensure that the first element is 0 (or the left truncation
#'        distance) and the last is the distance to the end of the furthest bin.
#'        (Default \code{NULL}, no binning.)
#'        Note that if \code{data} has columns \code{distbegin} and 
#'        \code{distend} then these will be used as bins if \code{cutpoints}
#'        is not specified. If both are specified, \code{cutpoints} has
#'        precedence.
#' @param monotonicity should the detection function be constrained for 
#'        monotonicity weakly ("weak"), strictly ("strict") or not at all 
#'        ("none" or \code{FALSE}). See Montonicity, below. (Default 
#'        \code{FALSE}).
#' @param region.table \code{data.frame} with two columns: 
#'        \tabular{ll}{ \code{Region.Label} \tab label for the region\cr
#'                     \code{Area} \tab area of the region\cr} 
#'        \code{region.table} has one row for each stratum. If there is no 
#'        stratification then \code{region.table} has one entry with \code{Area}
#'        corresponding to the total survey area.
#' @param sample.table \code{data.frame} mapping the regions to the samples (
#'        i.e. transects). There are three columns:
#'        \tabular{ll}{\code{Sample.Label} \tab label for the sample\cr
#'                     \code{Region.Label} \tab label for the region that the
#'                          sample belongs to.\cr
#'                     \code{Effort} \tab the effort expended in that sample 
#'                          (e.g. transect length).\cr} 
#' @param obs.table \code{data.frame} mapping the individual observations 
#'        (objects) to regions and samples. There should be three columns:
#'        \tabular{ll}{\code{object} \tab \cr
#'                     \code{Region.Label} \tab label for the region that the
#'                          sample belongs to.\cr
#'                     \code{Sample.Label} \tab label for the sample\cr}
#' @param convert.units conversion between units for abundance estimation, 
#'        see "Units", below. (Defaults to 1, implying all of the units are 
#'        "correct" already.)
#'        
#'
#' @return
#' TKTKTK need to think about this!!
#' @section Details:
#'
#'  If abundance estimates are required the \code{data.frame}s \code{region.table},
#'  \code{sample.table} and \code{obs.table} must be supplied.
#'
# THIS IS STOLEN FROM mrds, sorry Jeff!
#' @section Units:   
#'  In extrapolating to the entire survey region it is important that
#'  the unit measurements be consistent or converted for consistency.
#'  A conversion factor can be specified with the \code{convert.units}
#'  variable.  The values of \code{Area} in \code{region.table}, must be made 
#'  consistent with the units for \code{Effort} in \code{sample.table} and the 
#'  units of \code{distance} in the \code{data.frame} that was analyzed.  It is
#'  easiest if the units of \code{Area} are the square of the units of 
#'  \code{Effort} and then it is only necessary to convert the units of 
#'  \code{distance} to the units of \code{Effort}. For example, if \code{Effort}
#'   was entered in kilometers and \code{Area} in square kilometers and 
#'  \code{distance} in meters then using \code{convert.units=0.001} would 
#'  convert meters to kilometers, density would be expressed in square 
#'  kilometers which would then be consistent with units for \code{Area}.  
#'  However, they can all be in different units as long as the appropriate 
#'  composite value for \code{convert.units} is chosen.  Abundance for a survey 
#'  region can be expressed as: \code{A*N/a} where \code{A} is \code{Area} for 
#'  the survey region, \code{N} is the abundance in the covered (sampled) 
#'  region, and \code{a} is the area of the sampled region and is in units of 
#'  \code{Effort * distance}.  The sampled region \code{a} is multiplied by 
#'   \code{convert.units}, so it should be chosen such that the result is in 
#'  the same units as \code{Area}.  For example, if \code{Effort} was entered 
#'  in kilometers, \code{Area|} in hectares (100m x 100m) and \code{distance} 
#'  in meters, then using \code{convert.units=10} will convert \code{a} to 
#'  units of hectares (100 to convert meters to 100 meters for distance and 
#'  .1 to convert km to 100m units).
#'
#' @section Monotonicity:
#'   BLURB
#' @author David L. Miller
#' @export
#'
#' @examples
#'
#'  # An example from mrds, the golf tee data.
#'  library(distance)
#'  data(book.tee.data)
#'  tee.data<-book.tee.data$book.tee.dataframe
#'  ds.model<-ds(tee.data,4)
#'  summary(ds.model)
#'
#'  # same model, but calculating abundance
#'  # need to supply the region, sample and observation tables
#'  region<-book.tee.data$book.tee.region
#'  samples<-book.tee.data$book.tee.samples
#'  obs<-book.tee.data$book.tee.obs
#' 
#'  ds.dht.model<-ds(tee.data,4,region.table=region,
#'               sample.table=samples,obs.table=obs)
#'  summary(ds.dht.model)
#'
#'  # specify order 2 cosine adjustments
#'  ds.model.cos2<-ds(tee.data,4,adjustment="cos",order=2)
#'  summary(ds.model.cos2)
#'
#'  # specify order 2 and 4 cosine adjustments
#'  ds.model.cos24<-ds(tee.data,4,adjustment="cos",order=c(2,4))
#'  summary(ds.model.cos24)
#'
#'  # truncate the largest 10% of the data and fit only a hazard-rate
#'  # detection function
#'  ds.model.hr.trunc<-ds(tee.data,truncation="10%",key="hr",adjustment=NULL)
#'  summary(ds.model.hr.trunc)
#'
ds<-function(data, truncation=NULL, transect="line", formula=~1, key="hn",
             adjustment="cos", order=NULL, scale="scale", cutpoints=NULL,
             monotonicity=FALSE,
             region.table=NULL,sample.table=NULL,obs.table=NULL,
             convert.units=1){
  
  # this routine just creates a call to mrds, it's not very exciting
  # or fancy, it does do a lot of error checking though

  ## make sure that the data are in the right format first
  if(is.null(data$distance)){
    stop("Your data must have a column called 'distance'!")
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
  
  # truncation
  if(is.null(truncation)){
    stop("Please supply truncation distance or percentage.")
  }else{
    # if we have left truncation too...
    if(is.list(truncation)){
      if((any(names(truncation)=="left") & any(names(truncation)=="right")) &
          length(truncation)==2){

        # check for each of left and right that we have % or distance...
        # left
        if(is.double(truncation$left) & length(truncation$left)==1){
          left <- truncation$left
        }else if(is.character(truncation$left) & length(truncation$left)==1){
          # % string to number
          truncation$left<-as.numeric(sub("%","",truncation$left))
          left <- quantile(data$distance,probs=(truncation$left/100))
        }else{
          stop("Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
        }
        # right
        if(is.double(truncation$right) & length(truncation$right)==1){
          width <- truncation$right
        }else if(is.character(truncation$right) & length(truncation$right)==1){
          # % string to number
          truncation$right<-as.numeric(sub("%","",truncation$right))
          width <- quantile(data$distance,probs=1-(truncation$right/100))
        }else{
          stop("Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
        }
      }else{
        stop("Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
      }

    # just right truncation
    }else if(is.double(truncation) & length(truncation)==1){
      width <- truncation
      left <- NULL
    }else if(is.character(truncation) & length(truncation)==1){
      # % string to number
      truncation<-as.numeric(sub("%","",truncation))
      width <- quantile(data$distance,probs=1-(truncation/100))
      left <- NULL
    }else{
      stop("Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
    }
  }

  # transect type 
  if(transect=="line"){
    point <- FALSE
  }else if(transect=="point"){
    point <- TRUE
  }else{
    stop("Only \"point\" or \"line\" transects may be supplied.")
  }

  # key and adjustments
  if(!(key %in% c("hn","hr"))){
    stop("key function must be \"hn\" or \"hr\".\n")
  }
  if(!is.null(adjustment)){
    if(!(adjustment %in% c("cos","herm","poly"))){
      stop("adjustment terms must be one of NULL, \"cos\", \"herm\" or \"poly\".\n")
    }
  }
  if(!is.null(adjustment)){
    if(!is.null(order)){
      aic.search <- FALSE
      if(any(order != ceiling(order))){
          stop("Adjustment orders must be integers.")
      } 
      
      # check for each adjustment type
      order<-sort(order)
      if(adjustment=="herm" | adjustment=="poly"){
        if(any(order/2 != ceiling(order/2))){
          stop("Adjustment orders must be even for Hermite and simple polynomials.")
        }
      }
      if(adjustment=="herm" | adjustment=="cos"){
        if(any(order==1)){
          stop("Adjustment orders for Hermite polynomials and cosines must start at 2.")
        }
      }
    }else{
      aic.search <- TRUE
      max.order <- 3

      # this is according to p. 47 of IDS.
      if(adjustment=="poly"){
        order <- seq(1,max.order)
      }else{
        order <- seq(2,max.order)
      }
      if(adjustment=="herm" | adjustment=="poly"){
        order <- 2*order
        order <- order[order<=max.order]
      }
    }
  }else{
    aic.search<-FALSE
  }

  # binning
  if(is.null(cutpoints)){
    if(any(names(data)=="distend") & any(names(data)=="distbegin")){
      warning("No cutpoints specified but distbegin and distend are columns in data. Doing a binned analysis...")
      binned <- TRUE
    }else{
      binned <- FALSE
    }
  }else{
    # make sure that the first bin starts 0 or left
    if(!is.null(left)){
      if(cutpoints[1]!=left){
        stop("The first cutpoint must be 0 or the left truncation distance!")
      }
    }else if(cutpoints[1]!=0){
      stop("The first cutpoint must be 0 or the left truncation distance!")
    }

    # remove distbegin and distend if they already exist
    if(any(names(data)=="distend") & any(names(data)=="distbegin")){
      warning("data already has distend and distbegin columns, removing them and appling binning as specified by cutpoints.")
      data$distend<-NULL
      data$distbegin<-NULL
    }
    # send off to create.bins to make the correct columns in data
    data <- create.bins(data,cutpoints)
    binned <- TRUE
  }

  # monotonicity
  if(monotonicity=="none" | !monotonicity){
    mono<-FALSE
    mono.strict<-FALSE
  }else if(monotonicity=="weak"){
    mono<-TRUE
    mono.struct<-FALSE
  }else if(monotonicity=="strict"){
    mono<-TRUE
    mono.struct<-TRUE
  }else{
    stop("monotonicity must be one of \"none\", FALSE, \"weak\" or \"strict\".")
  } 
    
  ### Actually fit some models here
  
  models<-list()
  aics<-c()

  # construct the meta data object...
  meta.data <- list(width = width,point = point,binned = binned, 
                    mono=mono, mono.strict=mono.strict)
  if(!is.null(left)){
    meta.data$left<-left
  }

  # if we are doing an AIC-based search then, create the indices for the
  # for loop to work along, else just give the length of the order object
  if(aic.search){
    for.ind <- seq(along=order)
  }else if(!is.null(adjustment)){
    for.ind <- length(order)
  }else{
    for.ind <- 1
  }

  # loop over the orders of adjustments
  for(i in for.ind){ 
    # construct model formulae


    # CDS model
    if(as.character(formula)[2]=="1"){
      model.formula <- paste("~cds(key =\"", key,"\", formula = ~1",sep="")
    # MCDS model
    }else{
      model.formula <- paste("~mcds(key = \"",key,"\",", 
                                   "formula =~",as.character(formula)[2],sep="") 
    }

    # adjustments?
    if(!is.null(adjustment)){
      if(length(order[1:i])>1){
        order.str<-paste("c(",paste(order[1:i],collapse=","),")",sep="")
      }else{
        order.str<-order[1:i]
      }
      model.formula<-paste(model.formula,",", 
                           "adj.series=\"",adjustment,
                           "\",adj.order=",order.str,",",
                           "adj.scale=\"",scale,"\"",sep="")
    }

    model.formula<-paste(model.formula,")",sep="")


    # actually fit a model
    models[[i]]<-ddf(dsmodel = as.formula(model.formula),data = data, 
               method = "ds", meta.data = meta.data) 

    aics<-c(aics,models[[i]]$criterion)
  }

  # for AIC search, find the best fitting model
  if(aic.search){
    model<-models[[which.min(aics)]]
  }else{
    # otherwise just push the model back
    model<-models[[i]]
  }

  ## Now calculate abundance/density using dht()
  if(!is.null(region.table) & !is.null(sample.table) & !is.null(obs.table)){

    # from ?dht:
    # For animals observed in tight clusters, that estimator gives the 
    # abundance of groups (‘group=TRUE’ in ‘options’) and the abundance of 
    # individuals is estimated as s_1/p_1 + s_2/p_2 + ... + s_n/p_n, where 
    # s_i is the size (e.g., number of animals in the group) of each 
    # observation(‘group=FALSE’ in ‘options’).

    dht.res<-dht(model,region.table,sample.table,obs.table,
                 options=list(varflag=0,group=TRUE,
                              convert.units=convert.units),se=TRUE)
  }else{
    # if no information on the survey area was supplied just return 
    # the detection function stuff
    dht.res<-NULL

    cat("No survey area information supplied, only estimating detection function.\n")
  }

  # construct return object
  ret.obj<-list(dsmodel=model,dht=dht.res)

  # give it some class
  class(ret.obj)<-"distance"

  return(ret.obj)

}
