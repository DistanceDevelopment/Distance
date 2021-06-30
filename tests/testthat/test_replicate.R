# test replicate stratification option

context("Replicate")

manual_rep <- function(effort, density){

  tot.effort <- sum(effort)
  pool.over.grids <- sum(effort/tot.effort * density)

  var.pooled <- sum(effort*(density - pool.over.grids)^2) /
                  (tot.effort*(length(density) -1))
  se.pooled <- sqrt(var.pooled)
  cv.pooled <- se.pooled / pool.over.grids

  C <- exp(qnorm(1-0.025)*(sqrt(log(1+cv.pooled^2))))
  bh.90.LCI <- pool.over.grids/C
  bh.90.UCI <- pool.over.grids*C
  return(c(pool.over.grids, se.pooled, cv.pooled, bh.90.LCI, bh.90.UCI))
}

devtools::load_all()

test_that("Replicate works: minke", {
  skip_on_cran()

  data(minke)
  minke$Area <- 715316
  mink.cov <- ds(minke, key="hr", formula=~Region.Label)
  mink.dht2 <- dht2(mink.cov, flatfile = minke,
                    strat_formula=~Region.Label,
                    stratification = "replicate", total_area = 715316)
  den_res <- attr(mink.dht2, "density")

  # do manual calculations
  # density
  manual <- manual_rep(den_res$Effort[1:2], den_res$Density[1:2])
  # abundance
  manual2 <- manual_rep(mink.dht2$Effort[1:2], mink.dht2$Abundance[1:2])

  # test
  expect_equal(manual, as.numeric(den_res[3, c("Density", "Density_se",
                                               "Density_CV", "LCI", "UCI")]))
  expect_equal(manual2, as.numeric(mink.dht2[3, c("Abundance", "Abundance_se",
                                                  "Abundance_CV", "LCI",
                                                  "UCI")]))
})

test_that("Replicate works: Savannah sparrows", {
  skip_on_cran()

  data("Savannah_sparrow_1981")
  ss81 <- Savannah_sparrow_1981

  # relevel data for distance comparability
  ss81$Region.Label <- relevel(ss81$Region.Label, ref=4)


  cf <- convert_units("meter", NULL, "hectare")
  #it <- ds(ss81, key="hr", formula = ~Region.Label, convert.units = cf,
  #         transect="point", scale ="width",
  #         initial=list(scale=c(log(32.18), -0.1456, 0.2502, 0.8775E-01),
  #                      shape=log(4.549)))

  ss81d <- unflatten(ss81)$data
  # fit in ddf to fix initial values
  result <- ddf(dsmodel=~mcds(key="hr", formula=~Region.Label), data=ss81d,
                method="ds",
                meta.data=list(width=max(ss81d$distance), point=TRUE),
                control=list(initial=list(scale=c(log(32.18), -0.1456,
                                                  0.2502, 0.8775E-01),
                                          shape=log(4.549)),
                             nofit=TRUE))

  what <- dht2(result, flatfile = ss81, convert_units = cf,
               strat_formula=~Region.Label,
               stratification = "replicate")
  den_res <- attr(what, "density")
  manual <- manual_rep(den_res$Effort[1:4], den_res$Density[1:4])

  expect_equal(manual, as.numeric(den_res[5, c("Density", "Density_se",
                                               "Density_CV", "LCI", "UCI")]))
})
