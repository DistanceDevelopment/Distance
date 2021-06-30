# test replicate stratification option

context("Replicate")

test_that("Replicates: Savannah sparrows", {
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
  result <- ddf(dsmodel=~mcds(key="hr", formula=~1), data=ss81d,
                method="ds",
                meta.data=list(width=max(ss81d$distance), point=TRUE),
                control=list(initial=list(scale=log(33.23), shape=log(4.262)),
                             nofit=TRUE))

  # for density estimation
  ss81$Area <- 0

  what <- dht2(result, flatfile = ss81, convert_units = cf,
               strat_formula=~Region.Label,
               stratification = "replicate")

  expect_equal(what$Abundance, c(0.63466, 0.67093, 1.1424, 0.92479, 0.84319),
               tol=1e-3)

  expect_equal(what$Abundance_CV, c(18.88, 19.78, 16.67, 17.09, 16.88)/100,
               tol=1e-3)


  expect_equal(what$df[1:4], c(163.41, 156.26, 187.75, 182.20), tol=1e-1)
  expect_equal(what$df[5], 6.16, tol=1e-2)

  expect_equal(what$LCI, c(0.43861, 0.45567, 0.82406, 0.66168, 0.56098),
               tol=1e-3)

  expect_equal(what$UCI, c(0.91835, 0.98788, 1.5837, 1.2925, 1.2674),
               tol=1e-3)


# Distance results
#                         Estimate      %CV     df     95% Confidence Interval
#                        ------------------------------------------------------
# Stratum: PASTURE 0                                         
# Hazard/Cosine          
#                 D      0.63466       18.88   163.41 0.43861      0.91835    
#                 N          0.00000
# Stratum: PASTURE 1                                         
# Hazard/Cosine          
#                 D      0.67093       19.78   156.26 0.45567      0.98788    
#                 N          0.00000
# Stratum: PASTURE 2                                         
# Hazard/Cosine          
#                 D       1.1424       16.67   187.75 0.82406       1.5837    
#                 N          0.00000
# Stratum: PASTURE 3                                         
# Hazard/Cosine          
#                 D      0.92479       17.09   182.20 0.66168       1.2925    
#                 N          0.00000
#
# Pooled Estimates:
#                         Estimate      %CV     df     95% Confidence Interval
#                        ------------------------------------------------------
#                 D      0.84319       16.88     6.16 0.56098       1.2674    
#                 N          0.00000
#
})
