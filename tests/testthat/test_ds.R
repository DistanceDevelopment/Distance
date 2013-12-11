# test ds()

library(Distance)
library(mrds)

context("Testing ds()")


test_that("Input errors are thrown correctly",{

  data(book.tee.data)
  egdata<-book.tee.data$book.tee.dataframe

  # what if we don't supply any distance data?
  expect_error(ds())
  expect_error(ds(list()))

  # data but no truncation?
  expect_error(ds(egdata))

  # incorrect key definition?
  expect_error(ds(egdata,4,key="bananas"))

  # incorrect adjustment definition?
  expect_error(ds(egdata,4,key="hn",adjustment="bananas"))

  # first cutpoint not zero when no left truncation
  expect_error(ds(egdata,4,cutpoints=c(2,3,4)))

})


test_that("Simple models",{
  # data setup  
  data(book.tee.data)
  egdata<-book.tee.data$book.tee.dataframe

  test<-ds(egdata,4)
#'  # same model, but calculating abundance
#'  # need to supply the region, sample and observation tables
#'  region<-book.tee.data$book.tee.region
#'  samples<-book.tee.data$book.tee.samples
#'  obs<-book.tee.data$book.tee.obs
#' 
#'  ds.dht.model<-ds(tee.data,4,region.table=region,
#'               sample.table=samples,obs.table=obs)
#'  summary(ds.dht.model)
#'
#'  # specify order 2 cosine adjustments
#'  ds.model.cos2<-ds(tee.data,4,adjustment="cos",order=2)
#'  summary(ds.model.cos2)
#'
#'  # specify order 2 and 4 cosine adjustments
#'  ds.model.cos24<-ds(tee.data,4,adjustment="cos",order=c(2,4))
#'  summary(ds.model.cos24)
#'
#'  # truncate the largest 10% of the data and fit only a hazard-rate
#'  # detection function
#'  ds.model.hr.trunc<-ds(tee.data,truncation.percentage=10,key="hr",adjustment=NULL)
#'  summary(ds.model.hr.trunc)
  
})
