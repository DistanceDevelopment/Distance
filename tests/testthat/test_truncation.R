# test ds()

par.tol <- 1e-5
N.tol <- 1e-3
lnl.tol <- 1e-4

context("Testing truncation")

# data setup
data(book.tee.data)
egdata <- book.tee.data$book.tee.dataframe
egdata <- egdata[!duplicated(egdata$object),]

test_that("Truncation errors are thrown correctly",{

  # NULL truncation
  expect_error(ds(egdata, truncation=NULL),
                  "Please supply truncation distance or percentage.")

  # percentage truncation with bins
  expect_error(ds(egdata, truncation="10%",
                  cutpoints=seq(1, max(egdata$distance), length.out=10)),
                  "Truncation cannot be supplied as a percentage with binned data.")
  
  # wrong format - vector not list for 2 values
  expect_error(ds(egdata, truncation=c(left = 0, right = "10")),
               "Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\", each being a single value.")

  # too many left truncations
  expect_error(ds(egdata, truncation=list(left=c(0, 0), right=3.5)),
               "Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\", each being a single value.")
  
  # too many right truncations
  expect_error(ds(egdata, truncation=list(left=0, right=c(3.5, 3.5))),
               "Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\", each being a single value.")

  # too many characterleft truncations
  expect_error(ds(egdata, truncation=list(left=c("0%", "0%"), right=3.5)),
               "Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
  # too many character right truncations
  expect_error(ds(egdata, truncation=list(left=0, right=c("3.5%", "3.5%"))),
               "Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\", each being a single value.")

  # too many truncations
  expect_error(ds(egdata, truncation=list(left=0, right="3.5%", boop="3.5%")),
               "Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\", each being a single value.")

  # left and cutpoints[1] don't match
  expect_error(ds(egdata, truncation=list(left=0, right=max(egdata$distance)),
                  cutpoints=seq(1, max(egdata$distance), length.out=10)),
               "The first cutpoint must be 0 or the left truncation distance!")

  # message("data already has distend and distbegin columns, removing them and applying binning as specified by cutpoints.")
  
  # percentage truncation with no percent signs
  expect_warning(eg.1 <- ds(egdata, truncation=list(left = 0, right = "10")),
               "Truncation values supplied as characters will be interpreted as % truncation values.")
  
  eg.2 <- ds(egdata, truncation=list(left = 0, right = "10%"))
  
  # Should be the same
  expect_equal(eg.1$ddf, eg.2$ddf)
  # Check its calculated the quantile
  expect_equal(quantile(egdata$distance, probs = 0.9), eg.1$ddf$meta.data$width)
  expect_equal(quantile(egdata$distance, probs = 0.9), eg.2$ddf$meta.data$width)
  
  expect_warning(eg.3 <- ds(egdata, truncation=list(left = "10", right = "10%")),
                 "Truncation values supplied as characters will be interpreted as % truncation values.")
  # Check left truncation as a percentage
  expect_equal(quantile(egdata$distance, probs = 0.1), eg.3$ddf$meta.data$left)
  
  expect_warning(ds(egdata, truncation="10"),
                 "Truncation values supplied as characters will be interpreted as % truncation values.")
  
  # Check integer and double work the same
  eg.4 <- ds(egdata, truncation = 3L)
  eg.5 <- ds(egdata, truncation = 3)
  expect_equal(eg.4$ddf$fitted, eg.5$ddf$fitted)
  expect_equal(eg.4$ddf$Nhat, eg.5$ddf$Nhat)

})


