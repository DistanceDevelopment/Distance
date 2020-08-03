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
                  "Truncation cannot be supplied as a percentage with binned data")

  # too many left truncations
  expect_error(ds(egdata, truncation=list(left=c(0, 0), right=3.5)),
               "Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
  # too many right truncations
  expect_error(ds(egdata, truncation=list(left=0, right=c(3.5, 3.5))),
               "Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")

  # too many characterleft truncations
  expect_error(ds(egdata, truncation=list(left=c("0", "0"), right=3.5)),
               "Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
  # too many character right truncations
  expect_error(ds(egdata, truncation=list(left=0, right=c("3.5", "3.5"))),
               "Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")

  # too many truncations
  expect_error(ds(egdata, truncation=list(left=0, right="3.5", boop="3.5")),
               "Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")

  # left and cutpoints[1] don't match
  expect_error(ds(egdata, truncation=list(left=0, right=max(egdata$distance)),
                  cutpoints=seq(1, max(egdata$distance), length.out=10)),
               "The first cutpoint must be 0 or the left truncation distance!")

  # message("data already has distend and distbegin columns, removing them and appling binning as specified by cutpoints.")

})


