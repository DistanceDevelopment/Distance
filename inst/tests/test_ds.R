# test ds()

library(distance)
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
  expect_error(ds(egdata),4,key="bananas")

  # incorrect adjustment definition?
  expect_error(ds(egdata),4,key="hn",adjustment="bananas")

})


test_that("Simple models",{
  # data setup  
  data(book.tee.data)
  egdata<-book.tee.data$book.tee.dataframe

  test<-ds(egdata,4)
  
})
