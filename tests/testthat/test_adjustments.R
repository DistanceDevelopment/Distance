# test ds()

par.tol <- 1e-5
N.tol <- 1e-3
lnl.tol <- 1e-4

context("Adjustment terms")

# data setup
data(book.tee.data)
egdata <- book.tee.data$book.tee.dataframe
egdata <- egdata[!duplicated(egdata$object),]


test_that("adjustments expand correctly",{
  skip_on_cran()

  egdata <- egdata[egdata$observer==1,]

  # hn + cos(2)
  expect_equal(suppressWarnings(summary(ds(egdata, 4, key="hn",
                                           nadj=1))$ddf$name.message),
               "half-normal key function with cosine(2) adjustments")

  # hn + cos(2,3,4,5)
  expect_equal(suppressWarnings(summary(ds(egdata, 4, key="hn",
                                           nadj=4))$ddf$name.message),
               "half-normal key function with cosine(2,3,4,5) adjustments")

  #unif + cos(1,2)
  expect_equal(suppressWarnings(summary(ds(egdata, 4, key="unif",
                                           nadj=2))$ddf$name.message),
               "uniform key function with cosine(1,2) adjustments")

  #unif + poly(2)
  expect_equal(suppressWarnings(summary(ds(egdata, 4, key="unif",
                                           adjustment="poly",
                                           nadj=1))$ddf$name.message),
               "uniform key function with simple polynomial(2) adjustments")

  #unif + herm(2)
  expect_equal(suppressWarnings(summary(ds(egdata, 4, key="unif",
                                           adjustment="herm",
                                           nadj=1))$ddf$name.message),
               "uniform key function with Hermite(2) adjustments")

  # hn + cos(2,3)
  expect_equal(suppressWarnings(summary(ds(egdata, 4, key="hn",
                                           order=2:3))$ddf$name.message),
               "half-normal key function with cosine(2,3) adjustments")
})

test_that("adjustments orders start correctly",{
  skip_on_cran()

  # hn+poly starts at 4
  expect_message(suppressWarnings(ds(egdata, trunc=4, key="hn", adj="poly")),
                 "Fitting half-normal key function with simple polynomial\\(4\\) adjustments")
  # hn+cos starts at 2
  expect_message(suppressWarnings(ds(egdata, trunc=4, key="hn", adj="cos")),
                 "Fitting half-normal key function with cosine\\(2\\) adjustments")
  # hn+herm starts at 4
  expect_message(suppressWarnings(ds(egdata, trunc=4, key="hn", adj="herm")),
                 "Fitting half-normal key function with Hermite\\(4\\) adjustments")

  # hr+poly starts at 4
  expect_message(suppressWarnings(ds(egdata, trunc=4, key="hr", adj="poly")),
                 "Fitting hazard-rate key function with simple polynomial\\(4\\) adjustments")
  # hr+cos starts at 2
  expect_message(suppressWarnings(ds(egdata, trunc=4, key="hr", adj="cos")),
                 "Fitting hazard-rate key function with cosine\\(2\\) adjustments")
  # hr+herm starts at 4
  expect_message(suppressWarnings(ds(egdata, trunc=4, key="hr", adj="herm")),
                 "Fitting hazard-rate key function with Hermite\\(4\\) adjustments")

  # unif+poly starts at 2
  expect_message(suppressWarnings(ds(egdata, trunc=4, key="unif",
                                     adj="poly", max_adjustments=1)),
                 "Fitting uniform key function with simple polynomial\\(2\\) adjustments")
  # unif+cos starts at 1
  expect_message(suppressWarnings(ds(egdata, trunc=4, key="unif",
                                     adj="cos", max_adjustments=1)),
                 "Fitting uniform key function with cosine\\(1\\) adjustments")
  # unif+herm starts at 2
  expect_message(suppressWarnings(ds(egdata, trunc=4, key="unif",
                                     adj="herm", max_adjustments=1)),
                 "Fitting uniform key function with Hermite\\(2\\) adjustments")

})

# max adjustments arg
test_that("max.adjustments works",{
  skip_on_cran()

  egdata <- egdata[egdata$observer==1,]

  # setting max.adjustments=0 gives no adjustments
  expect_equal(summary(ds(egdata, 4, key="hn", max_adjustments=0,
                          adjustment="cos"))$ddf$name.message,
               "half-normal key function")

  # some delicious stake
  data(stake77)
  dists <- stake77$PD[stake77$Obs2==1]
  dists <- c(dists, dists[dists>10])
  dists <- c(dists, dists[dists<5])
  dists <- c(dists, dists[dists<5])

  # ignore warnings below from monotonicity checks, don't care about that here
  expect_equal(summary(suppressWarnings(
                ds(dists, 20, key="hn", max_adjustments=3,
                   adjustment="cos")))$ddf$name.message,
               "half-normal key function with cosine(2,3,4) adjustments")

  expect_equal(summary(suppressWarnings(
                ds(dists, 20, key="hn", max_adjustments=2,
                   adjustment="cos")))$ddf$name.message,
               "half-normal key function with cosine(2,3) adjustments")

  expect_equal(summary(suppressWarnings(
                ds(dists, 20, key="hn", max_adjustments=1,
                   adjustment="cos")))$ddf$name.message,
               "half-normal key function with cosine(2) adjustments")

  expect_equal(summary(suppressWarnings(
                ds(dists, 20, key="hn", max_adjustments=6,
                   adjustment="cos")))$ddf$name.message,
               "half-normal key function with cosine(2,3,4) adjustments")
})

test_that("errors thrown",{

  egdata <- egdata[egdata$observer==1,]

  # nadj and length(order) don't match
  expect_error(suppressWarnings(summary(ds(egdata, 4, key="hn", order=c(2,3),
                                           nadj=1))$ddf$name.message),
               "The number of adjustment orders specified in 'order' must match 'nadj'")

})
