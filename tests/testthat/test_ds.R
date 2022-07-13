# test ds()

par.tol <- 1e-5
N.tol <- 1e-3
lnl.tol <- 1e-4

context("Testing ds()")

# data setup
data(book.tee.data)
egdata <- book.tee.data$book.tee.dataframe
egdata <- egdata[!duplicated(egdata$object),]

test_that("Input errors are thrown correctly",{

  # what if we don't supply any distance data?
  expect_error(ds(),
               "argument \"data\" is missing, with no default")
  expect_error(ds(list()),
               "Your data must \\(at least\\) have a column called 'distance' or 'distbegin' and 'distend'!")

  # incorrect data types
  expect_error(ds("distance"),
               "data is not a data.frame nor are the values numeric")
  expect_error(ds(c(TRUE,FALSE)),
               "data is not a data.frame nor are the values numeric")

  # incorrect key definition?
  expect_error(ds(egdata,4,key="bananas"),
               suppressMessages(message("'arg' should be one of \"hn\", \"hr\", \"unif\"")))

  # incorrect adjustment definition?
  expect_error(ds(egdata,4,key="hn",adjustment="bananas"),
               suppressMessages(message("'arg' should be one of \"cos\", \"herm\", \"poly\"")))

  # uniform with covariates?
  expect_error(ds(egdata,4,key="unif",formula=~size),
               "Can't use uniform key with covariates.")

  # use a covariate that doesn't exist
  expect_error(ds(egdata,4,key="hn",formula=~banana),
               "Variable\\(s\\): banana are in the model formula but not in the data.")

  # observer is a reserved column
  egdata$observer <- c("Fred", "Larry")
  expect_error(ds(egdata,4,key="hn"),
               "'observer' is a reserved column name, please rename this column.")
  egdata$observer <- 1:3
  expect_error(ds(egdata,4,key="hn"),
               "'observer' is a reserved column name, please rename this column.")

})


test_that("binning works", {
  # first cutpoint not zero when no left truncation
  expect_error(ds(egdata,4,cutpoints=c(2,3,4)),
               "The first cutpoint must be 0 or the left truncation distance!")

  tst_distances <- data.frame(distance = c(0, 0, 0, 10, 50, 70, 110))
  expect_equal(as.vector(table(create_bins(tst_distances,
                                           c(0, 10, 65, 200))$distbegin)),
               c(4, 1, 2))

  # per bug #108
  data("wren_snapshot")

  binned <- expect_warning(create_bins(wren_snapshot,
                                       c(0, 10, 20, 30, 40, 60, 80, 100)),
                           "Some distances were outside bins and have been removed.")
  expect_equal(as.vector(table(binned$distbegin)), c(3, 9, 19, 46, 28, 11))
})

test_that("Simple models work",{
  skip_on_cran()

  #Test that ds can deal with a numeric vector as input
  distances <- c(1.02, 0.89, 0.21, 1.83, 0.09, 1.34, 2.1, 0.98, 1.8, 0.32,
                 0.83, 1.6, 0.92, 0.66, 0.31, 0.55, 1.13, 0.5, 0.46, 1)
  numeric.test <- suppressMessages(ds(distances))
  data.frame.test <- suppressMessages(ds(data.frame(distance = distances)))
  # Passing in distances as a vector should be identical to passing them in a
  # data.frame
  expect_equal(numeric.test$ddf$par, data.frame.test$ddf$par)

  # same model, but calculating abundance
  # need to supply the region, sample and observation tables
  region <- book.tee.data$book.tee.region
  samples <- book.tee.data$book.tee.samples
  obs <- book.tee.data$book.tee.obs

  egdata <- egdata[egdata$observer==1,]

  # half-normal key should get selected
  ds.dht.model <- suppressMessages(ds(egdata,4,region_table=region,
                                      sample_table=samples,obs_table=obs))
  # pars and lnl
  expect_equal(ds.dht.model$ddf$par, 0.6632435, tol=par.tol)
  expect_equal(ds.dht.model$ddf$lnl, -154.5692, tol=lnl.tol)
  expect_equal(ds.dht.model$dht$individuals$N$Estimate[3], 652.0909, tol=N.tol)

  # specify order 2 cosine adjustments
  ds.model.cos2<-suppressMessages(ds(egdata,4,adjustment="cos",order=2,
                                     region_table=region, sample_table=samples,
                                     obs_table=obs,monotonicity=FALSE))
  # pars and lnl
  #result <- ddf(dsmodel=~mcds(key="hn", formula=~1, adj.series="cos",
  #                            adj.order=2), data=egdata, method="ds",
  #              meta.data=list(width=4))
  tp <- c(0.66068510, -0.01592872)
  expect_equal(unname(ds.model.cos2$ddf$par), tp, tol=par.tol)
  expect_equal(ds.model.cos2$ddf$lnl, -154.5619, tol=1e-6)
  expect_equal(ds.model.cos2$dht$individuals$N$Estimate[3], 642.9442, tol=N.tol)

  # specify order 2 and 4 cosine adjustments
  ds.model.cos24<-suppressMessages(suppressWarnings(ds(egdata,4,
                     adjustment="cos",order=c(2,4),
                     region_table=region, sample_table=samples, obs_table=obs,
                     monotonicity=FALSE)))
  tp <- c(0.661391356, -0.008769883, -0.041465153)
  expect_equal(unname(ds.model.cos24$ddf$par), tp, tol=par.tol)
  expect_equal(ds.model.cos24$ddf$lnl, -154.5084, tol=lnl.tol)
  expect_equal(ds.model.cos24$dht$individuals$N$Estimate[3],620.0747,tol=N.tol)

  # hazard rate
  ds.model.hr<-suppressMessages(ds(egdata,4,key="hr",
                  adjustment=NULL, region_table=region,
                  sample_table=samples, obs_table=obs))
  #result <- ddf(dsmodel=~mcds(key="hr", formula=~1), data=egdata, method="ds",
  #              meta.data=list(width=4))
  tp <- c(0.9375891, 0.7245641)
  names(tp) <- c("p1","p2")
  expect_equal(ds.model.hr$ddf$par, tp, tol=par.tol)
  expect_equal(ds.model.hr$ddf$lnl, -155.2894, tol=lnl.tol)
  expect_equal(ds.model.hr$dht$individuals$N$Estimate[3], 591.8193, tol=N.tol)

})


test_that("Uniform does work after all",{
  skip_on_cran()

  egdata <- egdata[egdata$observer==1,]

  par.tol <- 1e-4
  # should select unif+cos(1)
  dd <- suppressMessages(ds(egdata,4,key="unif"))
  expect_equal(dd$ddf$par,  0.7384736, tol=par.tol)

  # try to fit with unif+cos(1,2)
  dd <- suppressMessages(ds(egdata,4,key="unif",order=c(1,2)))
  expect_equal(unname(dd$ddf$par), c(0.7050144, -0.1056291), tol=par.tol)

})


test_that("Truncation is handled",{
  skip_on_cran()

  egdata <- egdata[egdata$observer==1,]

  # setting the truncation is correct
  expect_equal(suppressMessages(ds(egdata, 4, key="hn",
                                   order=0))$ddf$meta.data$width, 4)

  # largest observed distance
  expect_equal(suppressMessages(ds(egdata, key="hn",
                                   order=0))$ddf$meta.data$width, 3.84)

  # remove observations after 3.8
  egdata <- egdata[egdata$distance <= 3.8, ]

  # largest cutpoint
  expect_equal(suppressMessages(ds(egdata,key="hn",order=0,
                             cutpoints=c(0,1,2,3,3.8)))$ddf$meta.data$width,3.8)

  # largest bin
  bin.data <- create_bins(egdata,c(0,1,2,3,3.8))
  expect_equal(suppressMessages(ds(bin.data, key="hn",
                                   order=0))$ddf$meta.data$width, 3.8)

})

# reported by Len Thomas 20 August
test_that("Percentage truncation works when distances are missing",{
  skip_on_cran()

  data(minke)

  expect_equal(ds(minke, truncation="15%", adjustment=NULL)$ddf$criterion,
               -8.1705496, tol=1e-5)
})


# just distend and distbegin can be supplied
test_that("just distend and distbegin can be supplied", {
  # make some data
  bin.data <- create_bins(egdata, c(0, 1, 2, 3, 4))
  bin.data$distance <- NULL
  expect_message(ds.model <- ds(bin.data, 4, monotonicity=FALSE, key="hn",
                                adjustment=NULL),
                 "^Columns \"distbegin\" and \"distend\" in data: performing a binned analysis....*")
})

test_that("cutpoints work with flatfile", {
  #https://github.com/DistanceDevelopment/Distance/issues/116
  skip_on_cran()
  distbegin <- c(rep(0,2), rep(1,36), rep(2, 29), rep(3, 17), rep(4,62),
                 rep(5,27),rep(6,11), rep(7,7), rep(NA,102))
  distend <- c(rep(1,2), rep(2,36), rep(3, 29), rep(4, 17), rep(5,62),
               rep(6,27),rep(7,11), rep(8,7), rep(NA,102))
  pigs <- data.frame(distance = (distbegin+distend)/2,
                     Region.Label = sample(c(rep("A", 200), rep("B", 93)),
                                           293, replace=TRUE),
                     Area         = 0)

  units <- convert_units("meter", NULL, "square kilometer")

  hn0_bushpig <- ds(pigs, transect = "point", key = "hn", adjustment = NULL,
                    convert_units = units,  cutpoints = c(seq(0,8,1)),
                    formula = ~Region.Label)

  hn0_bushpig <- ds(pigs, transect = "point", key = "hn", adjustment = NULL,
                    convert_units = units,  cutpoints = c(seq(0,8,1)),
                    formula = ~Region.Label)

  pigs2 <- data.frame(distance  = (distbegin+distend)/2,
                      distbegin = distbegin,
                      distend   = distend,
                      Region.Label = pigs$Region.Label,
                      Area         = 0)

  hn1_bushpig <- ds(pigs2, transect = "point", key = "hn", adjustment = NULL,
                    convert_units = units, formula = ~Region.Label)

  expect_equal(hn0_bushpig$ddf$par, hn1_bushpig$ddf$par, tol=par.tol)
  expect_equal(hn0_bushpig$ddf$lnl, hn1_bushpig$ddf$lnl, tol=lnl.tol)

})

# warnings when bad models get fitted
test_that("warnings of bad models get thrown", {
  skip_on_cran()

  # data with spike for which hazard-rate estimates
  # 0 shape par
  # based on data posted on distance list but rescaled and edited for anonymity
  dat <- c(0, 0, 0, 0, 0.008, 0.008, 0.013, 0.022, 0.027, 0.05, 0.058, 0.137,
           0.172, 0.2, 0.2, 0.2, 0.2, 0.25, 0.26, 0.29, 0.34, 0.36)
  expect_warning(ds(dat, key="hr", adjustment=NULL),
               "Estimated hazard-rate scale parameter close to 0 \\(on log scale\\). Possible problem in data \\(e.g., spike near zero distance\\).")

})
