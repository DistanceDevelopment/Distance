# variance exercise tests

context("variance")

test_that("variance 2",{

  data(Systematic_variance_2)

  # how to convert units
  convert.units <- Systematic_variance_2_units
  cu <- convert.units[, 3]
  cu <- 1/(cu[3]/(cu[1]*cu[2]))
  # first fit a model
  sysvar_df <- ds(Systematic_variance_2, adjustment="cos", convert_units=cu)

  systematic_var2 <- Systematic_variance_2
  unflat <- unflatten(systematic_var2)

  # make variance groupings
  vargroups <- data.frame(Sample.Label = 1:20,
                          grouping     = sort(rep(1:10, 2)))
  unflat$sample.table <- merge(unflat$sample.table, vargroups,
                               by="Sample.Label")

  # Remove these labels as they confuse the merge()
  unflat$sample.table$Region.Label <- NULL
  unflat$obs.table$Region.Label <- NULL
  unflat$sample.table$Sample.Label <-
    as.numeric(unflat$sample.table$Sample.Label)
  unflat$obs.table$Sample.Label <- as.numeric(unflat$obs.table$Sample.Label)
  unflat$data$Sample.Label <- as.numeric(unflat$data$Sample.Label)

  # fit the detection function
  sysvar_df <- ds(unflat$data, adjustment="cos", convert_units=cu)

  ## no stratification
  #Nhat_nostrat_eff <- dht2(sysvar_df, transects=unflat$sample.table,
  #                             geo_strat=unflat$region.table, observations=unflat$obs.table,
  #                             strat_formula=~1, convert_units=cu, er_est="R3")

  # using the Fewster et al 2009, "O2" estimator "Post-stratify, grouping adjacent pairs of samplers"
  expect_warning(Nhat_O2 <- dht2(sysvar_df, transects=unflat$sample.table,
                                 geo_strat=unflat$region.table,
                                 observations=unflat$obs.table,
                                 strat_formula=~1, convert_units=cu,
                                 er_est="O2"),
                 "Using the O2 encounter rate variance estimator, assuming that sorting on Sample.Label is meaningful")

  lr <- Nhat_O2[nrow(Nhat_O2), , drop=FALSE]
  expect_equal(lr$Abundance, 1022, tol=1e-1)
  expect_equal(lr$Abundance_CV, 8.64/100, tol=1e-2)
  expect_equal(lr$LCI, 858, tol=1e-1)
  expect_equal(lr$UCI, 1219, tol=1e-1)
  expect_equal(lr$df, 75.87, tol=1e-2)


# stratified estimate, weighting by effort
# (this is doing the effort weighted variance, combining the individual
# variance estimates)
# Nhat_strat_eff <- dht2(sysvar_df, transects=unflat$sample.table,
#                        geo_strat=unflat$region.table,
#                        observations=unflat$obs.table,
#                        strat_formula=~grouping, convert_units=cu,
#                        er_est="R3", stratification="effort_sum")
})
