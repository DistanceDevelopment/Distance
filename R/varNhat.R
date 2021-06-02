# this code is HORRIBLE
# Calculate the variance contribution of the detection function
#  to abundance estimates
varNhat <- function(data, model){

  # format the data
  # definitely a bad idea, relies on dplyr's internal representation
  # faff here is to detect if there is an NA grouping and drop it
  # this has all the group data in it
  grps <- attr(data, "groups")
  grps$.rows <- NULL
  # are there any NAs?
  grps <- apply(grps, 1, function(x) any(is.na(x)))
  # just get the indices that we want
  vardat_str <- attr(data, "groups")[!grps, , drop=FALSE]
  ind <- vardat_str$.rows

  # extract the area and covered area for each group
  area <- rep(NA, length(ind))
  covered_area <- rep(NA, length(ind))
  for(i in 1:length(ind)){
    idata <- data[ind[[i]], ]
    area[i] <- idata$Area[1]
# TODO: this is BAD
    covered_area[i] <- sum(idata$Covered_area[!duplicated(idata$Sample.Label)])

  }

  # function to calculate Nhat
  dhtd <- function(par, data, model, area, covered_area, ind){
    # result store
    res <- rep(NA, length(ind))

    # set par
    model$par <- par

    # calculate Nc
    data$p <- predict(model, newdata=as.data.frame(data), integrate=TRUE,
                      compute=TRUE)$fitted
    data$Nc <- (data$size/data$p)/data$rate

    # calculate Nhat per region
    for(i in 1:length(ind)){
      # need to na.rm to remove transects without obs
      res[i] <- (area[i]/covered_area[i]) * sum(data$Nc[ind[[i]]], na.rm=TRUE)
    }
    res
  }

  # get variance-covariance matrix for the detection function
  vcov <- solvecov(model$hessian)$inv

  # remove the rows where there were no observations
  data <- data[!is.na(data$object), ]

  # calculate variance
  dm <- DeltaMethod(model$par, dhtd, vcov, sqrt(.Machine$double.eps),
                    model=model, data=data, area=area,
                    covered_area=covered_area, ind=ind)
  attr(dm, "vardat_str") <- vardat_str

  ret <- list(Nhat=dm)

  attr(ret, "vardat_str") <- vardat_str
  return(ret)
}
