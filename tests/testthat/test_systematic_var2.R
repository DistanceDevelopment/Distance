# Systematic_variance_2

# test numerical tollerance
tol <- 1e-6

# make some results to check against
data(Systematic_variance_2)

Systematic_variance_2$Region.Label <- Systematic_variance_2$Region
dat <- unflatten(Systematic_variance_2)
obs.table <- dat$obs.table
region.table <- dat$region.table
sample.table <- dat$sample.table

context("systematic variance 2")


make_old_summ_cluster <- function(object){

  # print out as dht
  object$Region <- as.factor(object$Region)
  object$CoveredArea <- object$Covered_area
  object$se.ER <- sqrt(object$ER_var)
  object$cv.ER <- object$ER_CV
  object$se.mean <- sqrt(object$group_var)
  object$mean.size <- object$group_mean

  class(object) <- "data.frame"

  summ <- object[, c("Region", "Area", "CoveredArea", "Effort", "n", "k",
                     "ER", "se.ER", "cv.ER"), drop=FALSE]

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
  object[, c("Label", "Estimate", "se", "cv", "lcl", "ucl", "df"), drop=FALSE]
}



# ds fit
convert.units <- Systematic_variance_2_units
cu <- convert.units[, 3]
cu <- 1/(cu[3]/(cu[1]*cu[2]))


test_that("no strat",{

  df <- ds(Systematic_variance_2, convert_units=cu)


  # fiddle with region labels
  obs.table$Region.Label <- NULL
  sample.table$Region <- sample.table$Region.Label
  region.table$Region <- region.table$Region.Label
  sample.table$Region.Label <- NULL
  region.table$Region.Label <- NULL
  region.table$Region.Label <- NULL

  Systematic_variance_2$Region <- Systematic_variance_2$Region.Label
  Systematic_variance_2$Region.Label <- NULL

  # now do a fancy thing
  fs_st1 <- dht2(df$ddf, obs.table, sample.table, region.table,
                 strat_formula=~Region, convert_units=cu)

  # test
# work around stupid dht bug
df$dht$individuals$summary$k <- df$dht$individuals$summary$k[1]
  expect_equal(df$dht$individuals$summary,
               make_old_summ_cluster(fs_st1), tolerance=tol,
               check.attributes=FALSE)

  old <- df$dht$individuals$N
  old$Estimate <- as.numeric(old$Estimate)
  old$cv <- as.numeric(old$cv)
  old$ucl <- as.numeric(old$ucl)
  old$lcl <- as.numeric(old$lcl)
  old$df <- as.numeric(old$df)
  old$Label <- factor("Default")

  expect_equal(old, make_old_abund_individual(fs_st1), tolerance=tol)

})

