# cluster exercise

tol <- 1e-3

data(sikadeer)

cu <- sikadeer_units[, 3]
cu <- (cu[1]*cu[2])/cu[3]

context("sika deer")

#test_that("no truncation", {
# Effort        :    5.000000    
# # samples     :    37
# Width         :    200.0000    
# # observations:  1921
#
# Model  4
#    Half-normal key, k(y) = Exp(-y**2/(2*A(1)**2))
#    Cosine adjustments of order(s) :  2, 3, 4
#
#
#              Point        Standard    Percent Coef.        95% Percent
#  Parameter   Estimate       Error      of Variation     Confidence Interval
#  ---------  -----------  -----------  --------------  ----------------------
#    D         37.403       8.0972          21.65       24.355       57.441    
#    N         3400.0       736.04          21.65       2214.0       5221.0    
#  ---------  -----------  -----------  --------------  ----------------------
#
# Measurement Units                
# ---------------------------------
# Density: Numbers/Sq. kilometers 
#     ESW: centimeters    
#
# Component Percentages of Var(D)
# -------------------------------
# Detection probability   :   5.1
# Encounter rate          :  81.4
# Decay rate              :  13.6
#})

test_that("10% truncation", {

  dat <- unflatten(sikadeer)
  dat <- dat$data
  # fit the exact function fitted by Distance
  result <- ddf(dsmodel = ~mcds(key = "hn", formula = ~1,
                                adj.series="cos", adj.order=2),
                data = dat, method = "ds", meta.data = list(width = 175),
                control=list(initial=list(adjustment=0.08493070,
                                          scale=log(121.1144)),
                             nofit=TRUE))
  df <- result

  # multipliers
  mult <- list(creation = data.frame(rate = 25,
                                     SE   = 0),
               decay    = data.frame(rate = 163,
                                     SE   = 13))

  # weight by effort since we have repeats
  expect_warning(deer_ests <- dht2(df, flatfile=sikadeer,
                    strat_formula=~Region.Label,
                    convert_units=cu, multipliers=mult,
                    stratification="effort_sum",
                    total_area=13.9),
                 "One or more strata have only one transect, cannot calculate empirical encounter rate variance")

  deer_ests$ER_CV[is.na(deer_ests$ER_CV)] <- 0
  expect_equal(deer_ests$ER_CV[-nrow(deer_ests)],
               c(16.75, 24.13, 22.56, 46.95, 0, 49.61, 0, 0)/100, tol=tol)

  expect_equal(deer_ests$Abundance,
               c(1027.0, 383.00, 34.000, 29.000, 210.00, 126.00, 17.000,
                 69.000, 497), tol=1e-2)
  expect_equal(deer_ests$Abundance_CV,
               c(18.95, 25.71, 24.24, 47.78, 8.87, 50.39, 8.87, 8.87,
                 15.79)/100, tol=tol)
  expect_equal(deer_ests$df,
               c(19.18, 11.88, 2.74, 4.26, 47281.58, 2.07, 47281.58, 47281.58,
                 20.22), tol=tol)
  expect_equal(deer_ests$LCI,
               c(694, 220, 15, 9, 176, 18, 15, 58, 358), tol=1e-2)
  expect_equal(deer_ests$UCI,
               c(1521, 666, 76, 99, 249, 867, 21, 83, 689), tol=tol)
})
