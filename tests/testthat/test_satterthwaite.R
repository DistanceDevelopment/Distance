# test Satterthwaite degrees of freedom calculation


context("Satterthwaite degrees of freedom")


tol <- 1e-3

test_that("group size, no strat", {
  # add group sizes to minke
  groupsizes <- c(1, 1, 2, 4, 2, 6, 1, 3, 2, 5, 1, 2, 3, 2, 2, 1, 1, 2, 1, 2, 1,
                  1, 2, 2, NA, NA, 2, 2, 11, 1, 1, 2, NA, NA, NA, NA, 1, 1, 2, 2,
                  1, 2, 2, 2, 2, 1, 8, 10, 5, 7, 3, 1, 2, 2, 2, 10, 1, 1, 2, 2, 1,
                  2, NA, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 2, 3, 2, 1, 2, 1,
                  10, 1, 7, 2, 1, 1, 1, 1, 3, 1, NA, 1, 1, 1, 1, 1, NA)
  data(minke)
  minke$Region.Label <- "same"
  minke$Area <- 1000000
  minke$size <- groupsizes

  # fit and test dht
  simple <- ds(minke, key="hr", truncation=1.5, er_method=1, adjustment=NULL)

  expect_equal(simple$dht$individuals$D$Estimate, 0.57542E-01, tol=tol)
  expect_equal(simple$dht$clusters$D$Estimate, 0.25574E-01, tol=tol)
  expect_equal(simple$dht$individuals$D$df, 45.44, tol=tol)
  expect_equal(simple$dht$clusters$D$df, 34.49, tol=tol)
  expect_equal(simple$dht$individuals$D$cv, 27.89/100, tol=tol)
  expect_equal(simple$dht$clusters$D$cv, 25.98/100, tol=tol)


  # now with dht2
  d2 <- dht2(simple, flatfile=minke, strat_formula=~1, innes=FALSE)

  den_res <- attr(d2, "density")
  expect_equal(den_res$Density, 0.57542E-01, tol=tol)
  expect_equal(den_res$df, 45.44, tol=tol)
  expect_equal(den_res$Density_CV, 27.89/100, tol=tol)

  gr_den_res <- attr(attr(d2, "grouped"), "density")
  expect_equal(gr_den_res$Density, 0.25574E-01, tol=tol)
  expect_equal(gr_den_res$df, 34.49, tol=tol)
  expect_equal(gr_den_res$Density_CV, 25.98/100, tol=tol)

# distance results
# Effort        :    1842.790
# # samples     :    25
# Width         :    1.500000
# # observations:    88
#
# Model  1
#    Hazard Rate key, k(y) = 1 - Exp(-(y/A(1))**-A(2))
#
#
#              Point        Standard    Percent Coef.        95% Percent
#  Parameter   Estimate       Error      of Variation     Confidence Interval
#  ---------  -----------  -----------  --------------  ----------------------
#    DS       0.25574E-01  0.66433E-02      25.98      0.15219E-01  0.42975E-01
#    E(S)      2.2500      0.22872          10.17       1.8393       2.7524
#    D        0.57542E-01  0.16051E-01      27.89      0.33159E-01  0.99853E-01
#    N         41161.       11482.          27.89       23719.       71427.
#  ---------  -----------  -----------  --------------  ----------------------
#
# Measurement Units
# ---------------------------------
# Density: Numbers/Sq. kilometers
#     ESW: kilometers
#
# Component Percentages of Var(D)
# -------------------------------
# Detection probability   :  14.8
# Encounter rate          :  71.9
# Cluster size            :  13.3


#                         Estimate      %CV     df     95% Confidence Interval
#                        ------------------------------------------------------
# Hazard/Cosine          
#                 DS     0.25574E-01   25.98    34.49 0.15219E-01  0.42975E-01
#                 D      0.57542E-01   27.89    45.44 0.33159E-01  0.99853E-01
#                 N       41161.       27.89    45.44  23719.       71427.    

})




test_that("stratification works", {
  groupsizes <- c(1, 1, 2, 4, 2, 6, 1, 3, 2, 5, 1, 2, 3, 2, 2, 1, 1, 2, 1, 2,
                  1, 1, 2, 2, NA, NA, 2, 2, 11, 1, 1, 2, NA, NA, NA, NA, 1, 1, 2,
                  2, 1, 2, 2, 2, 2, 1, 8, 10, 5, 7, 3, 1, 2, 2, 2, 10, 1, 1, 2,
                  2, 1, 2, NA, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 2, 3, 2, 1,
                  2, 1, 10, 1, 7, 2, 1, 1, 1, 1, 3, 1, NA, 1, 1, 1, 1, 1, NA)
  data(minke)
  minke$size <- groupsizes

  simple2 <- ds(minke, key="hr", truncation=1.5, er_method=1, adjustment=NULL)

  expect_equal(simple2$dht$individuals$D$Estimate,
               c(0.44945E-01, 0.92867E-01, 0.50621E-01), tol=tol)
  expect_equal(simple2$dht$clusters$D$Estimate,
               c(0.19318E-01, 0.43117E-01, 0.22137E-01), tol=tol)
  expect_equal(simple2$dht$individuals$D$df, c(17.01, 28.87, 19.64), tol=1e-2)
  expect_equal(simple2$dht$clusters$D$df, c(12.97, 17.96, 15.26), tol=tol)
  expect_equal(simple2$dht$individuals$D$cv, c(40.82, 28.33, 33.13)/100, tol=tol)
  expect_equal(simple2$dht$clusters$D$cv, c(38.08, 24.91, 30.53)/100, tol=tol)



  # now with dht2
  d2 <- dht2(simple2, flatfile=minke, strat_formula=~Region.Label, innes=FALSE)
  den_res <- attr(d2, "density")
  expect_equal(den_res$Density,
               c(0.44945E-01, 0.92867E-01, 0.50621E-01), tol=tol)
  expect_equal(den_res$df, c(17.01, 28.87, 19.64), tol=1e-2)
  expect_equal(den_res$Density_CV, c(40.82, 28.33, 33.13)/100, tol=1e-2)

  gr_den_res <- attr(attr(d2, "grouped"), "density")
  expect_equal(gr_den_res$Density,
               c(0.19318E-01, 0.43117E-01, 0.22137E-01), tol=tol)
  expect_equal(gr_den_res$df, c(12.97, 17.96, 15.26), tol=tol)
  expect_equal(gr_den_res$Density_CV, c(38.08, 24.91, 30.53)/100, tol=1e-2)

#                         Estimate      %CV     df     95% Confidence Interval
#                        ------------------------------------------------------
# Stratum: North                                             
# Hazard/Cosine          
#                 DS     0.19318E-01   38.08    12.97 0.87222E-02  0.42787E-01
#                 D      0.44945E-01   40.82    17.01 0.19635E-01  0.10288    
#                 N       28341.       40.82    17.01  12381.       64874.    
# Stratum: South                                             
# Hazard/Cosine          
#                 DS     0.43117E-01   24.91    17.96 0.25747E-01  0.72205E-01
#                 D      0.92867E-01   28.33    28.87 0.52604E-01  0.16394    
#                 N       7869.0       28.33    28.87  4457.0       13892.    
#
#
# Pooled Estimates:
#                         Estimate      %CV     df     95% Confidence Interval
#                        ------------------------------------------------------
#                 DS     0.22137E-01   30.53    15.26 0.11727E-01  0.41789E-01
#                 D      0.50621E-01   33.13    19.64 0.25800E-01  0.99321E-01
#                 N       36210.       33.13    19.64  18455.       71046.    
})
