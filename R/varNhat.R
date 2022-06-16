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

  grp_dat <- data %>%
    select(Covered_area, Area, Sample.Label) %>%
    distinct() %>%
    summarize(Covered_area = sum(Covered_area),
              Area         = Area) %>%
    distinct()

  data$Covered_area <- NULL
  data$Area <- NULL

  data <- left_join(data, grp_dat)

  # function to calculate Nhat
  dhtd <- function(par, data, model){#, area, covered_area, ind){
    # set par
    model$par <- par

    ## calculate Nc
    data$p <- predict(model, newdata=as.data.frame(data), integrate=TRUE,
                      compute=TRUE)$fitted

    res <- data %>%
      mutate(Nc = (.data$size/.data$p)/.data$rate) %>%
      summarize(N = (.data$Area/.data$Covered_area) *
                     sum(.data$Nc, na.rm=TRUE)) %>%
      distinct()

    res$N
  }

  # get variance-covariance matrix for the detection function
  vcov <- solvecov(model$hessian)$inv

  # remove the rows where there were no observations
  data <- data[!is.na(data$object), ]

  # calculate variance
  dm <- DeltaMethod(model$par, dhtd, vcov, sqrt(.Machine$double.eps),
                    model=model, data=data)
  attr(dm, "vardat_str") <- vardat_str

  ret <- list(Nhat=dm)

  attr(ret, "vardat_str") <- vardat_str
  return(ret)
}
