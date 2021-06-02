# cue counting example

# test numerical tollerance
tol <- 1e-6

# make some results to check against
data(golftees)

# data fiddling
gtees <- golftees[golftees$observer==1 & golftees$detected==1, ]
gtees$size <- 1
dat <- unflatten(gtees)

trunc <- 4

make_old_summ_cluster <- function(object){

  # print out as dht
  object$Region <- as.factor(object$Region)
  object$CoveredArea <- object$Covered_area
  object$se.ER <- sqrt(object$ER_var)
  object$cv.ER <- object$ER_CV
  object$se.mean <- sqrt(object$group_var)
  object$mean.size <- object$group_mean

  class(object) <- "data.frame"

  summ <- object[, c("Region", "Area", "CoveredArea", "Effort", "n",
                     "ER", "se.ER", "cv.ER", "mean.size", "se.mean")]
  summ$se.ER[is.na(summ$se.ER) | is.nan(summ$se.ER)] <- 0
  summ$cv.ER[is.na(summ$cv.ER) | is.nan(summ$cv.ER)] <- 0
  return(summ)
}
make_old_abund_individual <- function(object){

  object$Label <- as.factor(object$Region)
  object$Region <- NULL
  object$Estimate <- object$Abundance
  object$cv <- object$Abundance_CV
  object$se <- object$Abundance_se
  object$lcl <- object$LCI
  object$ucl <- object$UCI
  object$df <- object$df
  class(object) <- "data.frame"
  object[, c("Label", "Estimate", "se", "cv", "lcl", "ucl", "df")]
}


context("golftees")

test_that("ER variance",{

  # output using old dht
  df <- ds(gtees, truncation=trunc, key="hn", adjustment=NULL)

  # now do a fancy thing
  dat$obs.table <- dat$obs.table[dat$obs.table$object %in% gtees$object, ]
  fs_st1 <- expect_warning(dht2(df$ddf, dat$obs.table, dat$sample.table,
                                dat$region.table, strat_formula=~Region.Label,
                                innes=FALSE),
                           "One or more strata have only one transect, cannot calculate empirical encounter rate variance")

  # tests
  oldres <- df$dht$clusters$summary
  oldres$k <- NULL
  oldres$mean.size <- df$dht$individuals$summary$mean.size
  oldres$se.mean <- df$dht$individuals$summary$se.mean
  expect_equal(oldres,
               make_old_summ_cluster(fs_st1), tolerance=tol)

  fs_st1$Region <- "Total"
  # TODO: this is very stupid and probably a bug in mrds
#  aa <- df$dht$individuals$N
#  aa[, -1] <- as.numeric(aa[, -1])
#  expect_equal(aa, make_old_abund_individual(fs_st1),
#               tolerance=tol, check.attributes=FALSE)

})
