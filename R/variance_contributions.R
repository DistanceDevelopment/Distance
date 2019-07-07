# what are the variance contributions for each stratum?
variance_contributions <- function(res){

  # extract the columns we want
  CV_cont <- data.frame(ER          = res$ER_CV,
                        Groups      = sqrt(res$group_var)/res$group_mean,
                        Multipliers = res$rate_CV,
                        Detection   = sqrt(res$p_var)/res$p_average)

  # remove Multipliers if not there
  if(all(is.na(CV_cont$Multipliers) | is.nan(CV_cont$Multipliers) |
         CV_cont$Multipliers == 0)) CV_cont$Multipliers <- NULL

  # remove group size if not there
  if(all(is.na(CV_cont$Groups) | is.nan(CV_cont$Groups) |
         CV_cont$Groups == 0)) CV_cont$Groups <- NULL

  # get the total
  CV_cont$Total <- sqrt(rowSums(CV_cont^2))

  # make that into percentages
  CV_cont <- (CV_cont^2/CV_cont[["Total"]]^2)*100
  CV_cont[["Total"]] <- NULL

  # zero ER contributions if only one sample
  CV_cont$ER[res$k==1] <- 0

  # sort and name
  CV_cont <- cbind(res[,1], CV_cont[order(names(CV_cont))])
  names(CV_cont)[1] <- names(res)[1]

  return(CV_cont)
}

