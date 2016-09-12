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
               "Your data must \\(at least\\) have a column called 'distance'!")
  
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

  # first cutpoint not zero when no left truncation
  expect_error(ds(egdata,4,cutpoints=c(2,3,4)),
               "The first cutpoint must be 0 or the left truncation distance!")

  # uniform with covariates?
  expect_error(ds(egdata,4,key="unif",formula=~size),
               "Can't use uniform key with covariates.")

  # use a covariate that doesn't exist
  expect_error(ds(egdata,4,key="hn",formula=~banana),
               "Variable\\(s\\): banana are in the model formula but not in the data.")


})


test_that("Simple models work",{
  
  #Test that ds can deal with a numeric vector as input
  distances <- c(1.02,0.89,0.21,1.83,0.09,1.34,2.1,0.98,1.8,0.32,0.83,1.6,0.92,0.66,0.31,0.55,1.13,0.5,0.46,1)
  numeric.test <- suppressMessages(ds(distances))
  data.frame.test <- suppressMessages(ds(data.frame(distance = distances)))
  #Passing in distances as a vector should be identical to passing them in a data.frame
  expect_equal(numeric.test$ddf$par, data.frame.test$ddf$par)

  # same model, but calculating abundance
  # need to supply the region, sample and observation tables
  region<-book.tee.data$book.tee.region
  samples<-book.tee.data$book.tee.samples
  obs<-book.tee.data$book.tee.obs

  egdata <- egdata[egdata$observer==1,]

  # half-normal key should get selected
  ds.dht.model <- suppressMessages(ds(egdata,4,region.table=region,
               sample.table=samples,obs.table=obs))
  # pars and lnl
  expect_equal(ds.dht.model$ddf$par, 0.6632435,tol=par.tol)
  expect_equal(ds.dht.model$ddf$lnl, -154.5692, tol=lnl.tol)
  expect_equal(ds.dht.model$dht$individuals$N$Estimate[3], 652.0909, tol=N.tol)

  # specify order 2 cosine adjustments
  ds.model.cos2<-suppressMessages(ds(egdata,4,adjustment="cos",order=2,
                                     region.table=region, sample.table=samples,
                                     obs.table=obs,monotonicity=FALSE))
  # pars and lnl
  #result <- ddf(dsmodel=~mcds(key="hn", formula=~1, adj.series="cos",
  #                            adj.order=2), data=egdata, method="ds",
  #              meta.data=list(width=4))
  tp <- c(0.66068510, -0.01592872)
  names(tp) <- c("X.Intercept.","V2")
  expect_equal(ds.model.cos2$ddf$par, tp, tol=par.tol)
  expect_equal(ds.model.cos2$ddf$lnl, -154.5619, tol=1e-6)
  expect_equal(ds.model.cos2$dht$individuals$N$Estimate[3], 642.9442, tol=N.tol)

  # specify order 2 and 4 cosine adjustments
  ds.model.cos24<-suppressMessages(ds(egdata,4,adjustment="cos",order=c(2,4),
                     region.table=region, sample.table=samples, obs.table=obs,
                     monotonicity=FALSE))
  tp <- c(0.661391356, -0.008769883, -0.041465153)
  names(tp) <- c("X.Intercept.","V2","V3")
  expect_equal(ds.model.cos24$ddf$par, tp, tol=par.tol)
  expect_equal(ds.model.cos24$ddf$lnl, -154.5084, tol=lnl.tol)
  expect_equal(ds.model.cos24$dht$individuals$N$Estimate[3],620.0747,tol=N.tol)

  # hazard rate
  ds.model.hr<-suppressMessages(ds(egdata,4,key="hr",
                  adjustment=NULL, region.table=region,
                  sample.table=samples, obs.table=obs))
  #result <- ddf(dsmodel=~mcds(key="hr", formula=~1), data=egdata, method="ds",
  #              meta.data=list(width=4))
  tp <- c(0.9375891, 0.7245641)
  names(tp) <- c("p1","p2")
  expect_equal(ds.model.hr$ddf$par, tp, tol=par.tol)
  expect_equal(ds.model.hr$ddf$lnl, -155.2894, tol=lnl.tol)
  expect_equal(ds.model.hr$dht$individuals$N$Estimate[3], 591.8193, tol=N.tol)

})


test_that("Uniform does work after all",{

  egdata <- egdata[egdata$observer==1,]

  par.tol <- 1e-4
  # should select unif+cos(1)
  dd <- suppressMessages(ds(egdata,4,key="unif"))
  expect_equal(dd$ddf$par,  0.7384736, tol=par.tol)

  # try to fit with unif+cos(1,2)
  dd <- suppressMessages(ds(egdata,4,key="unif",order=c(1,2)))
  expect_equal(dd$ddf$par, c(0.7050144, -0.1056291), tol=par.tol)


})



test_that("Truncation is handled",{

  egdata <- egdata[egdata$observer==1,]

  # setting the truncation is correct
  expect_equal(suppressMessages(ds(egdata,4,key="hn",order=0))$ddf$meta.data$width,4)

  # largest observed distance
  expect_equal(suppressMessages(ds(egdata,key="hn",order=0))$ddf$meta.data$width,3.84)


  # remove observations after 3.8
  egdata <- egdata[egdata$distance <= 3.8,]

  # largest cutpoint
  expect_equal(suppressMessages(ds(egdata,key="hn",order=0,
                  cutpoints=c(0,1,2,3,3.8)))$ddf$meta.data$width,3.8)

  # largest bin
  bin.data <- Distance:::create.bins(egdata,c(0,1,2,3,3.8))
  expect_equal(suppressMessages(ds(bin.data,key="hn",order=0))$ddf$meta.data$width,3.8)

})

test_that("adjustments expand correctly",{

  egdata <- egdata[egdata$observer==1,]

  # hn + cos(2)
  expect_equal(summary(ds(egdata,4,key="hn",order=2))$ddf$name.message,
               "half-normal key function with cosine(2) adjustments")

  # hn + cos(2,3)
  expect_equal(summary(ds(egdata,4,key="hn",order=2:3))$ddf$name.message,
               "half-normal key function with cosine(2,3) adjustments")

  # hn + cos(2,3,4,5)
  expect_equal(summary(ds(egdata,4,key="hn",order=5))$ddf$name.message,
               "half-normal key function with cosine(2,3,4,5) adjustments")

  #unif + cos(1,2)
  expect_equal(summary(ds(egdata,4,key="unif",order=2))$ddf$name.message,
               "uniform key function with cosine(1,2) adjustments")

})

# reported by Len Thomas 20 August
test_that("Percentage truncation works when distances are missing",{

  data(minke)

  expect_equal(ds(minke, truncation="15%", adjustment=NULL)$ddf$criterion,
               -8.1705496, tol=1e-5)
})
