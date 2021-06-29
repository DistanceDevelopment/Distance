# minke

# test numerical tollerance
tol <- 1e-6

# make some results to check against
data(minke, package="Distance")

minke$Region.Label <- as.factor(minke$Region.Label)

# set truncation
whale.trunc <- 1.5


# save no object # version
minke_noobj <- minke

# data fiddling
minke$object <- NA
minke$object[!is.na(minke$distance)] <- 1:nrow(minke[!is.na(minke$distance),])

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
                     "ER", "se.ER", "cv.ER")]#, "mean.size", "se.mean")]
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


context("minke")

test_that("South works",{
  skip_on_cran()

  whales <- minke[minke$Region.Label=="South",]
  dat <- unflatten(whales)

  # fiddle with region labels
  dat$obs$Region.Label <- NULL
  dat$sample$Region <- dat$sample$Region.Label
  dat$region$Region <- dat$region$Region.Label
  dat$sample$Region.Label <- NULL
  dat$region$Region.Label <- NULL

  # output using old dht
  whale.full.strat1 <- ds(whales, truncation=whale.trunc,
                          key="hr",  adjustment=NULL)
  whale.full.strat1$ddf$data$Region.Label <- NULL

  # now do a fancy thing
  fs_st1 <- with(dat, dht2(whale.full.strat1$ddf, obs, sample, region,
                     strat_formula=~Region))

  # test
  whale.full.strat1$dht$individuals$summary$k <-
    whale.full.strat1$dht$individuals$summary$k[1]
  expect_equal(whale.full.strat1$dht$individuals$summary,
               make_old_summ_cluster(fs_st1), tolerance=tol)

  old <- whale.full.strat1$dht$individuals$N
  old$Estimate <- as.numeric(old$Estimate)
  old$cv <- as.numeric(old$cv)
  old$ucl <- as.numeric(old$ucl)
  old$lcl <- as.numeric(old$lcl)
  old$df <- as.numeric(old$df)
  old$Label <- as.factor("South")
  expect_equal(old, make_old_abund_individual(fs_st1), tolerance=tol)

})

test_that("North works",{
  skip_on_cran()

  whales <- minke[minke$Region.Label=="North",]
  dat <- unflatten(whales)

  # fiddle with region labels
  dat$obs$Region.Label <- NULL
  dat$sample$Region <- dat$sample$Region.Label
  dat$region$Region <- dat$region$Region.Label
  dat$sample$Region.Label <- NULL
  dat$region$Region.Label <- NULL

  # output using old dht
  whale.full.strat1 <- ds(whales, truncation=whale.trunc,
                          key="hr",  adjustment=NULL)
  whale.full.strat1$ddf$data$Region.Label <- NULL

  # now do a fancy thing
  fs_st1 <- with(dat, dht2(whale.full.strat1$ddf, obs, sample, region,
                     strat_formula=~Region))

  # test
  whale.full.strat1$dht$individuals$summary$k <-
    whale.full.strat1$dht$individuals$summary$k[1]
  expect_equal(whale.full.strat1$dht$individuals$summary,
               make_old_summ_cluster(fs_st1), tolerance=tol)

  #TODO: fix this!
  old <- whale.full.strat1$dht$individuals$N
  old$Estimate <- as.numeric(old$Estimate)
  old$cv <- as.numeric(old$cv)
  old$ucl <- as.numeric(old$ucl)
  old$lcl <- as.numeric(old$lcl)
  old$df <- as.numeric(old$df)
  old$Label <- as.factor("North")
  expect_equal(old, make_old_abund_individual(fs_st1), tolerance=tol)

})

test_that("Combined works",{
  skip_on_cran()

  whales <- minke
  whales$Region.Label <- as.factor("Total")
  whales$Area <- sum(unique(whales$Area))
  dat <- unflatten(whales)

  # fiddle with region labels
  dat$obs$Region.Label <- NULL
  dat$sample$Region <- dat$sample$Region.Label
  dat$region$Region <- dat$region$Region.Label
  dat$sample$Region.Label <- NULL
  dat$region$Region.Label <- NULL

  # output using old dht
  whale.full.strat1 <- ds(whales, truncation=whale.trunc,
                          key="hr",  adjustment=NULL)
  whale.full.strat1$ddf$data$Region.Label <- NULL

  # now do a fancy thing
  fs_st1 <- with(dat, dht2(whale.full.strat1$ddf, obs, sample,
                               region, strat_formula=~Region))

  # test
  whale.full.strat1$dht$individuals$summary$k <-
    whale.full.strat1$dht$individuals$summary$k[1]
  expect_equal(whale.full.strat1$dht$individuals$summary,
               make_old_summ_cluster(fs_st1), tolerance=tol)

  #TODO: fix this!
  old <- whale.full.strat1$dht$individuals$N
  old$Estimate <- as.numeric(old$Estimate)
  old$cv <- as.numeric(old$cv)
  old$ucl <- as.numeric(old$ucl)
  old$lcl <- as.numeric(old$lcl)
  old$df <- as.numeric(old$df)
  old$Label <- as.factor(old$Label)
  expect_equal(old, make_old_abund_individual(fs_st1), tolerance=tol)

})


test_that("stratification", {
  skip_on_cran()

  whales <- minke
  # cutpoints
  cp <- seq(0, 1.5, length=8)

  expect_warning(wmr <- ds(whales, truncation=whale.trunc,
                           key="hr",  adjustment=NULL, cutpoints=cp),
                 "Some distances were outside bins and have been removed.")

  expect_warning(whale.full.strat1 <- ds(whales, truncation=whale.trunc,
                                         key="hr", adjustment=NULL,
                                         cutpoints=cp),
                 "Some distances were outside bins and have been removed.")

  result <- ddf(dsmodel = ~mcds(key = "hr", formula = ~1),
                data = whale.full.strat1$ddf$data, method = "ds",
                meta.data = list(width = 1.5, binned=TRUE, breaks=cp),
                control=list(initial=list(shape=log(2.483188),
                                          scale=log(0.7067555)),
                             nofit=TRUE))

  whale.full.strat1 <- result

  # now do a fancy thing
  fs_st1 <- dht2(whale.full.strat1, flatfile=whales,
                 strat_formula=~Region.Label)
#                         Estimate      %CV     df     95% Confidence Interval
#                        ------------------------------------------------------
# Stratum: 1. South                                          
# Hazard/Cosine          
#                 D      0.44471E-01   25.27    18.97 0.26419E-01  0.74858E-01
#                 N       3768.0       25.27    18.97  2239.0       6343.0    
# Stratum: 2. North                                          
# Hazard/Cosine          
#                 D      0.19925E-01   38.31    13.29 0.89727E-02  0.44247E-01
#                 N       12564.       38.31    13.29  5658.0       27901.    
#	Estimation Summary - Density&Abundance       	
#
# Pooled Estimates:
#                         Estimate      %CV     df     95% Confidence Interval
#                        ------------------------------------------------------
#                 D      0.22833E-01   30.82    15.83 0.12050E-01  0.43266E-01
#                 N       16333.       30.82    15.83  8619.0       30949.    

#                       Estimate      %CV     df     95% Confidence Interval
#                        ------------------------------------------------------
# Stratum: 1. South                                          
#                 n       39.000    
#                 k       13.000    
#                 L       484.41    
#                 n/L    0.80510E-01   22.48    12.00 0.49630E-01  0.13061    
#                 Left    0.0000
#                 Width   1.5000    
# Stratum: 2. North                                          
#                 n       49.000    
#                 k       12.000    
#                 L       1358.4    
#                 n/L    0.36072E-01   36.54    11.00 0.16551E-01  0.78620E-01
#                 Left    0.0000
#                 Width   1.5000    

expect_equal(fs_st1$Abundance, c(12564, 3768, 16333), tol=1e-3)
expect_equal(fs_st1$LCI, c(5658, 2239, 8619), tol=1e-3)
expect_equal(fs_st1$UCI, c(27901, 6343, 30949), tol=1e-3)
expect_equal(fs_st1$ER_CV[1:2], c(36.54, 22.48)/100, tol=1e-2)
expect_equal(fs_st1$Abundance_CV, c(38.31, 25.27, 30.82)/100, tol=1e-2)
expect_equal(fs_st1$df, c(13.29, 18.97, 15.83), tol=1e-2)

})


test_that("flat minke works", {
  skip_on_cran()

  # compare dht and dht2
  strat.spec <- ds(minke_noobj, key="hr", formula=~Region.Label, tru=1.5)

  fromdht2 <- dht2(strat.spec, flatfile=minke_noobj,
                   strat_formula=~Region.Label)

  expect_equal(strat.spec$dht$individuals$N$Estimate,
               fromdht2$Abundance)

  expect_equal(strat.spec$dht$individuals$N$se,
               fromdht2$Abundance_se, tol=1e-4)

})

test_that("minke dht vs dht2", {
  skip_on_cran()

  uf <- unflatten(minke)

  # compare dht and dht2
  mfit <- ds(minke, key="hr", formula=~Region.Label, tru=1.5)
  dht_old <- dht(mfit$ddf, uf$region.table, uf$sample.table,
                 options=list(ervar="R2", varflag=1))

  fromdht2 <- dht2(mfit, flatfile=minke,
                   strat_formula=~Region.Label)

  expect_equal(dht_old$individuals$N$Estimate,
               fromdht2$Abundance)
  #expect_equal(dht_old$individuals$summary$se.ER,
  #             sqrt(fromdht2$ER_var), tol=1e-4)

  expect_equal(dht_old$individuals$N$se,
               fromdht2$Abundance_se, tol=1e-4)


  # now with innes
  dht_old <- dht(mfit$ddf, uf$region.table, uf$sample.table,
                 options=list(ervar="R2", varflag=2))
  fromdht2 <- dht2(mfit, flatfile=minke, strat_formula=~Region.Label, innes=TRUE)

  expect_equal(dht_old$individuals$N$Estimate,
               fromdht2$Abundance)
  # note this will always be wrong because dht reports ER without Innes
  # but dht2 reports with Innes
  #expect_equal(dht_old$individuals$summary$se.ER,
  #             sqrt(fromdht2$ER_var), tol=1e-4)
  expect_equal(dht_old$individuals$N$se,
               fromdht2$Abundance_se, tol=1e-4)
})

test_that("area=0, flat minke works", {
  skip_on_cran()

  minke_noobj$Area <- 0
  # compare dht and dht2
  strat.spec <- ds(minke_noobj, key="hr", formula=~Region.Label, tru=1.5)

  fromdht2 <- dht2(strat.spec, flatfile=minke_noobj,
                   strat_formula=~Region.Label)

  expect_equal(strat.spec$dht$individuals$D$Estimate,
               fromdht2$Abundance)

  expect_equal(strat.spec$dht$individuals$D$se,
               fromdht2$Abundance_se, tol=1e-4)

  expect_equal(strat.spec$dht$individuals$summary$se.ER, sqrt(fromdht2$ER_var),
               tol=1e-3)


})


