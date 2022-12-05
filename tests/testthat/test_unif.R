library(Distance)
library(testthat)

context("Uniform models")

test_that("models with a uniform key function are handled correctly",{
  
  data(book.tee.data)
  tee.data <- subset(book.tee.data$book.tee.dataframe, observer==1)
  
  # Uniform with no adjustments Pa should be 1
  ds.model <- ds(tee.data, 4,  key = "unif", nadj = 0)
  #expect_true(all(ds.model$ddf$fitted == 1))
  expect_null(ds.model$ddf$par)
  
  ds.model <- ds(tee.data, 4,  key = "unif", nadj = 0)
  #expect_true(all(ds.model$ddf$fitted == 1))
  expect_null(ds.model$ddf$par)
  
  # Uniform with adjustments specified
  
  # Uniform with no adjustments and left truncation 
  
  
})