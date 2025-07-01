par.tol <- 1e-5
N.tol <- 1e-3
lnl.tol <- 1e-4

# Create test dataset
bin.data <- data.frame(Region.Label = "Region",
                       Area = 100,
                       Sample.Label = c(rep("A",10), rep("B",10), rep("C",13)),
                       Effort = 10,
                       distance = c(rep(5,11), rep(15,10), rep(35,9), rep(75,3)))


context("Binned Analyses")

test_that("binned analysis arguments are interpreted correctly",{
  skip_on_cran()
  
  # Analysis using distance and cutpoints
  bin.fit <- ds(bin.data,
                truncation = 100,
                formula = ~1,
                key = "hr",
                nadj = 0,
                cutpoints = c(0,10,20,50,100))
  
  #add distbegin and distend to the data
  BEbin.data <- create_bins(bin.data, cutpoints = c(0,10,20,50,100))
  BEbin.data$distance <- NULL
  
  # Should produce the same analysis
  bin.fit2 <- ds(BEbin.data,
                 truncation = 100,
                 formula = ~1,
                 key = "hr",
                 nadj = 0)
  
  library(testthat)
  expect_equal(bin.fit$ddf$par, bin.fit2$ddf$par)
  expect_equal(bin.fit$ddf$meta.data, bin.fit2$ddf$meta.data)
  
  # Test warning if user provides a distance column as well as distbegin and distend
  both.bin.data <- create_bins(bin.data, cutpoints = c(0,10,20,50,100))
  expect_warning(bin.fit3 <- ds(both.bin.data,
                                truncation = 100,
                                formula = ~1,
                                key = "hr",
                                nadj = 0),
                 "You have supplied both a 'distance' column and 'distbegin' and 'distend' columns in your data, the distance column will be removed and not used in these analyses.")
  
  expect_equal(bin.fit$ddf$meta.data, bin.fit3$ddf$meta.data)
  
  # Check there is an error if both distbegin and distend are not supplied together
  begin.bin.data <- create_bins(bin.data, cutpoints = c(0,10,20,50,100))
  begin.bin.data$distend <- NULL
  expect_error(ds(begin.bin.data, truncation = 100,
                  formula = ~1, key = "hr", nadj = 0),
               "You have provided either a 'distbegin' or 'distend' column in your dataset but not both. Please provide both or remove these and provide a distance column and use the cutpoint argument.")
  
  # Check that is it not using cutpoints when distbegin and distend are in the data
  expect_warning(bin.fit4 <- ds(BEbin.data, truncation = 100,
                                formula = ~1, key = "hr", nadj = 0, 
                                cutpoints = c(0,20,50,100)),
                 "Data has distend and distbegin columns, cutpoints argument will be ignored.")
  
  expect_equal(bin.fit$ddf$meta.data, bin.fit4$ddf$meta.data)
  
})
