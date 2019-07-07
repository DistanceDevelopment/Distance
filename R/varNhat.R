# this code is HORRIBLE
# Calculate the variance contribution of the detection function
#  to abundance estimates
varNhat <- function(data, model, strata){

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

  # remove data where we don't have observations, since this doesn't effect
  # the predictions...
  #data <- data[!is.na(data$object), ]

  # function to calculate Nhat
  dhtd <- function(par, data, model, area, covered_area, ind){
    # result store
    res <- rep(NA, length(ind))

    # set par
    model$par <- par

    # calculate Nc
    data$p <- predict(model, newdata=as.data.frame(data), integrate=TRUE, compute=TRUE)$fitted
    data$Nc <- (data$size/data$p)/data$rate

    # calculate Nhat per region
    for(i in 1:length(ind)){
      # need to na.rm to remove transects without obs
      res[i] <- (area[i]/covered_area[i]) * sum(data$Nc[ind[[i]]], na.rm=TRUE)
    }
    res
  }

#  # function to calculate P_a
#  phtd <- function(par, data, model, ind){
#    # result store
#    res <- rep(NA, length(ind))
#
#    # set par
#    model$par <- par
#
#    # calculate Nc
#    data$p <- predict(model, newdata=data, integrate=TRUE, compute=TRUE)$fitted
#    data$Nc <- 1/data$p
#
#    # calculate Nhat per region
#    for(i in 1:length(ind)){
#      # p = n/Nhat
#      res[i] <- sum(data$n_observations[ind[[i]]], na.rm=TRUE)/
#                 sum(data$Nc[ind[[i]]], na.rm=TRUE)
#    }
#    res
#  }


  # get variance-covariance matrix for the detection function
  vcov <- solvecov(model$hessian)$inv

  # calculate variance
  dm <- DeltaMethod(model$par, dhtd, vcov, sqrt(.Machine$double.eps),
                    model=model, data=data, area=area,
                    covered_area=covered_area, ind=ind)
  attr(dm, "vardat_str") <- vardat_str
#  # calculate variance
#  pm <- DeltaMethod(model$par, phtd, vcov, sqrt(.Machine$double.eps),
#                    model=model, data=data, ind=ind)
#
#  covar <- t(dm$partial) %*% vcov %*% pm$partial + pm$covar
#  average.p <- phtd(model$par, data, model, ind)
#  Nhat <- dhtd(model$par, data, model, area=area, covered_area=covered_area, ind)
#  var.pbar <- average.p^2*((sqrt(dm$variance)/Nhat)^2 + pm$var/data$n_observations^2-
#                                  2*covar/(data$n_observations*Nhat))
#

  ret <- list(Nhat=dm)#, phat=var.pbar)

  attr(ret, "vardat_str") <- vardat_str
  return(ret)
}
