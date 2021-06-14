#' Print abundance estimates
#' @export
#' @param x object of class `dht_result`
#' @param groups should abundance/density of groups be produced?
#' @param report should `"abundance"`, `"density"` or `"both"` be reported?
#' @param \dots unused
print.dht_result <- function(x, report="abundance", groups=FALSE, ...){

  # general information
  resulttype <- ifelse(report=="both", "abundance and density", report)
  substr(resulttype, 1, 1) <- toupper(substr(resulttype, 1, 1))
  cat(resulttype, "estimates from distance sampling\n")
  cat("Stratification :", attr(x, "stratification"), "\n")
  cat("Variance       :", paste0(attr(x, "ER_var")[1], ","),
      ifelse(attr(x, "ER_var")[3], "binomial",
             ifelse(attr(x, "ER_var")[2], "N/L", "n/L")), "\n")
  cat("Multipliers    :", ifelse(is.null(attr(x, "multipliers")),
                                 "none",
                                 paste(attr(x, "multipliers"), collapse=", ")),
      "\n")
  cat("Sample fraction :" , ifelse(attr(x, "sample_fraction")>1,
                                  "multiple", attr(x, "sample_fraction")), "\n")

  cat("\n\n")
  if(groups & !is.null(attr(x, "grouped"))){
    cat("Groups:\n\n")
    print(attr(x, "grouped"), report=report)
    cat("\n\n")
    cat("Individuals:\n\n")
  }

  object <- x
  # print out as dht
  object$CoveredArea <- object$Covered_area
  object$se.ER <- sqrt(object$ER_var)
  object$cv.ER <- sqrt(object$ER_var)/object$ER
  object$cv.ER[object$ER==0] <- 0

  # get the stratification labels
  stratum_labels <- attr(object, "stratum_labels")

  class(object) <- "data.frame"

  # print summary statistics
  cat("Summary statistics:\n")
  summ <- object[, c(stratum_labels, "Area", "CoveredArea", "Effort", "n",
                     "k", "ER", "se.ER", "cv.ER")]

  if(!is.null(attr(object, "summary"))){
    ss <- attr(object, "summary")
    ss$CoveredArea <- ss$Covered_area
    ss$se.ER <- sqrt(ss$ER_var)
    ss$cv.ER <- ss$ER_CV
    summ <- rbind(summ, ss[, c(stratum_labels, "Area", "CoveredArea",
                               "Effort", "n","k", "ER", "se.ER", "cv.ER")])
  }
  summ[,c("ER", "se.ER", "cv.ER")] <- round(summ[,c("ER", "se.ER", "cv.ER")], 3)
  print(summ, row.names=FALSE)
  cat("\n")

  # which columns?
  if(all(object$group_mean == 1 | is.na(object$group_mean) |
     is.nan(object$group_mean))){
    cols <- c(stratum_labels, "Estimate", "se", "cv", "LCI",
                   "UCI", "df")
  }else{
    object$group_se <- sqrt(object$group_var)
    cols <- c(stratum_labels, "Estimate", "se", "cv", "LCI",
                   "UCI", "df", "group_mean", "group_se")
  }
  if(report=="abundance" | report=="both" | attr(object, "density_only")){
    # print abundance table
    if(attr(object, "density_only")){
      round <- 4
      cat("Density estimates:\n")
    }else{
      round <- 0
      cat("Abundance estimates:\n")
    }
    object$Estimate <- object$Abundance
    object$cv <- object$Abundance_CV
    object$se <- object$Abundance_se
    ab <- object[, cols]

    ab[,c("Estimate", "LCI", "UCI")] <- round(ab[,c("Estimate", "LCI", "UCI")],
                                              round)
    ab[,c("se", "cv", "df")] <- round(ab[,c("se", "cv", "df")], 3)
    print(ab, row.names=FALSE)
    cat("\n")
  }
  # density estimates if requested
  if((report=="density" | report=="both") & !attr(object, "density_only")){
      round <- 4
    cat("Density estimates:\n")
    dobject <- attr(object, "density")

    dobject$Estimate <- dobject$Density
    dobject$cv <- dobject$Density_CV
    dobject$se <- dobject$Density_se
    dobject$group_se <- sqrt(dobject$group_var)
    ab <- dobject[, cols]

    ab[,c("Estimate", "LCI", "UCI")] <- round(ab[,c("Estimate", "LCI", "UCI")],
                                              round)
    ab[,c("se", "cv", "df")] <- round(ab[,c("se", "cv", "df")], 3)
    print(ab, row.names=FALSE)
    cat("\n")
  }
  cat("Component percentages of variance:\n")
  var_cont <- attr(object, "prop_var")
  var_cont[, -1] <- round(var_cont[, -1], 2)
  print(var_cont, row.names=FALSE)
  cat("\n")

  if(report=="abundance" & attr(object, "density_only")){
    cat("Can't report abundance, only density was estimated.\n")
  }

  invisible()
}
