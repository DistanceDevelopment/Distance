# test replicate stratification option

context("Replicate")

manual_rep <- function(effort, density){

  tot.effort <- sum(effort)
  pool.over.grids <- sum(effort/tot.effort * density)

  var.pooled <- sum((density - pool.over.grids)^2) / (length(density) -1)
  se.pooled <- sqrt(var.pooled)
  cv.pooled <- se.pooled / pool.over.grids

  C <- exp(1.96*(sqrt(log(1+cv.pooled^2))))
  bh.90.LCI <- pool.over.grids/C
  bh.90.UCI <- pool.over.grids*C
  return(c(pool.over.grids, se.pooled, cv.pooled, bh.90.LCI, bh.90.UCI))
}


test_that("Replicate works: minke", {
  skip_on_cran()

devtools::load_all()
  data(minke)
  mink.cov <- ds(minke, key="hr", formula=~Region.Label)
  mink.dht2 <- dht2(mink.cov, flatfile = minke,
                    strat_formula=~Region.Label,
                    stratification = "replicate", total_area = 715316)
  den_res <- attr(mink.dht2, "density")
  manual <- manual_rep(den_res$Effort[1:2], den_res$Density[1:2])

  expect_equal(manual, as.numeric(den_res[3, c("Density", "Density_se", "Density_CV", "LCI", "UCI")]))
})

test_that("Replicate works: Savannah sparrows", {
  skip_on_cran()

  data("Savannah_sparrow_1981")
  ss81 <- Savannah_sparrow_1981
  cf <- convert_units("meter", NULL, "hectare")
  it <- ds(ss81, key="hr", formula = ~Region.Label, convert.units = cf)
  what <- dht2(it, flatfile = ss81, convert_units = cf,
               strat_formula=~Region.Label,
               stratification = "replicate")
  den_res <- attr(what, "density")
  manual <- manual_rep(den_res$Effort[1:4], den_res$Density[1:4])

  expect_equal(manual, as.numeric(den_res[5, c("Density", "Density_se", "Density_CV", "LCI", "UCI")]))
})
