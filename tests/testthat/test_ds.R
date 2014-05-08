# test ds()

library(Distance)
library(mrds)

context("Testing ds()")

# data setup
data(book.tee.data)
egdata<-book.tee.data$book.tee.dataframe

test_that("Input errors are thrown correctly",{

  # what if we don't supply any distance data?
  expect_error(ds())
  expect_error(ds(list()))

  # incorrect key definition?
  expect_error(ds(egdata,4,key="bananas"))

  # incorrect adjustment definition?
  expect_error(ds(egdata,4,key="hn",adjustment="bananas"))

  # first cutpoint not zero when no left truncation
  expect_error(ds(egdata,4,cutpoints=c(2,3,4)))

  # uniform with covariates?
  expect_error(ds(egdata,4,key="unif",formula=~size))

})


test_that("Simple models work",{

  # same model, but calculating abundance
  # need to supply the region, sample and observation tables
  region<-book.tee.data$book.tee.region
  samples<-book.tee.data$book.tee.samples
  obs<-book.tee.data$book.tee.obs

  # half-normal key should get selected
  ds.dht.model<-ds(egdata,4,region.table=region,
               sample.table=samples,obs.table=obs)
  # pars and lnl
  expect_equal(ds.dht.model$ddf$par, 0.9229255)
  expect_equal(ds.dht.model$ddf$lnl, -215.1466, tol=1e-6)
  expect_equal(ds.dht.model$dht$individuals$N$Estimate[3], 712.6061, tol=1e-6)

  # specify order 2 cosine adjustments
  ds.model.cos2<-ds(egdata,4,adjustment="cos",order=2, region.table=region,
                    sample.table=samples,obs.table=obs,monotonicity=FALSE)
  # pars and lnl
  #result <- ddf(dsmodel=~mcds(key="hn", formula=~1, adj.series="cos",
  #                            adj.order=2), data=egdata, method="ds",
  #              meta.data=list(width=4))
  tp <- c(0.92098796, -0.04225286)
  names(tp) <- c("X.Intercept.","V2")
  expect_equal(ds.model.cos2$ddf$par, tp)
  expect_equal(ds.model.cos2$ddf$lnl, -215.0763, tol=1e-6)
  expect_equal(ds.model.cos2$dht$individuals$N$Estimate[3], 682.5572, tol=1e-6)

  # specify order 2 and 4 cosine adjustments
  ds.model.cos24<-ds(egdata,4,adjustment="cos",order=c(2,4),
                     region.table=region, sample.table=samples, obs.table=obs,
                     monotonicity=FALSE)
  tp <- c(0.92121582, -0.03712634, -0.03495348)
  names(tp) <- c("X.Intercept.","V2","V3")
  expect_equal(ds.model.cos24$ddf$par, tp)
  expect_equal(ds.model.cos24$ddf$lnl, -215.0277, tol=1e-6)
  expect_equal(ds.model.cos24$dht$individuals$N$Estimate[3], 661.1458, tol=1e-6)

  # hazard rate
  ds.model.hr<-ds(egdata,4,key="hr",
                  adjustment=NULL, region.table=region,
                  sample.table=samples, obs.table=obs)
  #result <- ddf(dsmodel=~mcds(key="hr", formula=~1), data=egdata, method="ds",
  #              meta.data=list(width=4))
  tp <- c(0.7833921, 0.9495860)
  names(tp) <- c("p1","p2")
  expect_equal(ds.model.hr$ddf$par, tp, tol=1e-7)
  expect_equal(ds.model.hr$ddf$lnl, -215.2384, tol=1e-6)
  expect_equal(ds.model.hr$dht$individuals$N$Estimate[3], 661.9913, tol=1e-6)

})


test_that("Truncation is handled",{

  # setting the truncation is correct
  expect_equal(ds(egdata,4,key="hn",order=0)$ddf$meta.data$width,4)

  # largest observed distance
  expect_equal(ds(egdata,key="hn",order=0)$ddf$meta.data$width,3.84)

  # largest cutpoint
  expect_equal(ds(egdata,key="hn",order=0,cutpoints=c(0,1,2,3,3.8))$ddf$meta.data$width,3.8)

  # largest bin
  bin.data <- Distance:::create.bins(egdata,c(0,1,2,3,3.8))
  expect_equal(ds(bin.data,key="hn",order=0)$ddf$meta.data$width,3.8)

})
