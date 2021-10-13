context("Testing bootdht")

test_that("area=0 with default summary",{

  skip_on_cran()

  data(amakihi)
  # subset for testing speed
  amakihi <- amakihi[1:300,]
  conv <- convert_units("meter", NULL, "hectare")
  surv <- ds(amakihi, transect="point", key="hr",
             formula=~Region.Label, convert_units=conv,
             truncation = 82.5, er_var="P3")
  expect_error(bootdht(surv, flatfile=amakihi, nboot=5),
               "No Area in flatfile, densities will be returned and the default summary function records only abundances. You need to write your own summary_fun.")

  amakihi$Area <- NULL
  expect_error(bootdht(surv, flatfile=amakihi, nboot=5),
               "No Area in flatfile, densities will be returned and the default summary function records only abundances. You need to write your own summary_fun.")

})


# generate some data to test on
dat <- data.frame(Sample.Label = 1:10)
dat$Region.Label <- c(rep(1, 5), rep(2, 5))
set.seed(123)
dat2 <- dat[sample(1:nrow(dat), nrow(dat), replace=TRUE), ]
dat2$object <- 1:nrow(dat2)
dat <- dat[!(dat$Sample.Label %in% dat2$Sample.Label),]
dat$object <- NA
dat <- rbind(dat, dat2)
dat$save.ind <- 1:nrow(dat)


test_that("resamples work - strata", {

  set.seed(123)
  obs <- bootdht_resample_data(dat, c("Region.Label"))

  expect_equal(obs$save.ind, c(1, 4, 5, 7, 9, 10, 2, 3, 6, 13, 8, 11, 12))
})


test_that("resamples work - sample", {

  set.seed(123)
  obs <- bootdht_resample_data(dat, c("Sample.Label"))

  expect_equal(obs$save.ind, c(7, 7, 4, 5, 4, 5, 7, 12, 8, 11, 2, 3, 6, 13))
})

test_that("resamples work - object", {

  set.seed(123)
  obs <- bootdht_resample_data(dat, c("object"))

  expect_equal(obs$save.ind, c(1, 4, 5, 6, 5, 2, 3, 6, 11, 9))
})

test_that("resamples work - strata/object", {

  set.seed(123)
  obs <- bootdht_resample_data(dat, c("Region.Label", "object"))

  expect_equal(obs$save.ind, c(1, 4, 5, 5, 5, 9, 2, 3, 13, 8, 11, 8))
})


test_that("resamples work - sample/object", {

  set.seed(123)
  obs <- bootdht_resample_data(dat, c("Sample.Label", "object"))

  expect_equal(obs$save.ind, c(7, 4, 4, 4, 5, 4, 12, 8, 2, 3, 6, 13))
})

test_that("resamples work - strata/sample/object", {

  set.seed(1223)
  obs <- bootdht_resample_data(dat, c("Region.Label", "Sample.Label", "object"))

  expect_equal(obs$save.ind, c(6, 1, 9, 6, 6, 2, 3, 8, 11, 6, 13, 8))
})
