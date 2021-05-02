# input checking

data(minke, package="Distance")
minke$object <- 1:nrow(minke)
dat <- unflatten(minke)
df <- ds(dat$data, truncation=1.5)
obs <- dat$obs
transects <- dat$sample.table
strata <- dat$region.table
strata$Stratum <- strata$Region.Label
strata$Region.Label <- NULL

context("input checking")


test_that("observation columns", {

  bad_obs <- obs
  colnames(bad_obs) <- c("boop", "blorp")
  expect_error(dht2(df, bad_obs, transects, strat_formula=~1),
               "^observations data must at least contain columns named")

  bad_obs <- obs
  colnames(bad_obs) <- c("object", "blorp")
  expect_error(dht2(df, bad_obs, transects, strat_formula=~1),
               "^observations data must at least contain columns named")

  bad_obs <- obs
  colnames(bad_obs) <- c("Sample.Label", "blorp")
  expect_error(dht2(df, bad_obs, transects, strat_formula=~1),
               "^observations data must at least contain columns named")

})


test_that("transect columns", {

  bad_transects <- transects
  colnames(bad_transects) <- c("boop", "blorp", "jumblethorp")
  expect_error(dht2(df, obs, bad_transects, strat_formula=~1),
               "^transects data must at least contain columns named")

  bad_transects <- transects
  colnames(bad_transects) <- c("Effort", "blorp", "jumblethorp")
  expect_error(dht2(df, obs, bad_transects, strat_formula=~1),
               "^transects data must at least contain columns named")

  bad_transects <- transects
  colnames(bad_transects) <- c("Sample.Label", "blorp", "jumblethorp")
  expect_error(dht2(df, obs, bad_transects, strat_formula=~1),
               "^transects data must at least contain columns named")

})


test_that("geo_strat columns", {

  bad_geo_strat <- strata
  colnames(bad_geo_strat) <- c("Stratum", "blorp")
  expect_error(dht2(df, obs, transects, bad_geo_strat, strat_formula=~Stratum),
               "^geo_strat data must at least contain columns named")

  bad_geo_strat <- strata
  colnames(bad_geo_strat) <- c("Area", "blorp")
  bad_transects <- transects
  colnames(bad_transects) <- c("Effort", "Sample.Label", "jumblethorp")
  # also need to drop that stratum from the transects
  expect_error(dht2(df, obs, bad_transects, bad_geo_strat,
                        strat_formula=~Stratum),
               "Stratification variable\\(s\\) \"Stratum\" not found in the data")

  # what if there are dupes?
  bad_geo_strat <- rbind(strata, strata)
  expect_error(dht2(df, obs, transects, bad_geo_strat,
                        strat_formula=~Stratum),
               "Inconsistent Areas/stratum labels in `geo_strat`")



})


test_that("stratification variable is in the data", {
  expect_error(dht2(df, obs, transects, strata, strat_formula=~foo),
               "Stratification variable(s) \"foo\" not found in the data",
               fixed=TRUE)
})


