test_that("mutlipliers errors",{
  data("minke")
  easy <- ds(data=minke, key="hr", truncation=1.5)
  # bad labels
  mult <- list(cration = data.frame(rate=c(0.41, 0.47),
                                     SE=c(0.07, 0.06)))
  expect_error(res_dens <- dht2(easy, flatfile=minke,
                                strat_formula = ~Region.Label,
                                stratification="geographical",
                                multipliers = mult),
               "^Multipliers must be named \"creation\" and \"decay\"")

  # rows without labels
  mult <- list(creation = data.frame(rate=c(0.41, 0.47),
                                     SE=c(0.07, 0.06)))
  expect_error(res_dens <- dht2(easy, flatfile=minke,
                                strat_formula = ~Region.Label,
                                stratification="geographical",
                                multipliers = mult),
               "^Multirow multipliers need column to link to the data")

  # labels but they aren't in the data
  mult <- list(creation = data.frame(rate=c(0.41, 0.47),
                                     SE=c(0.07, 0.06),
                                     Region.Labelz=c("North", "South")))
  expect_error(res_dens <- dht2(easy, flatfile=minke,
                                strat_formula = ~Region.Label,
                                stratification="geographical",
                                multipliers = mult),
               "^Multirow multipliers need column to link to the data")

})
