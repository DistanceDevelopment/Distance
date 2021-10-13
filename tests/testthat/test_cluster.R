# cluster tests

context("cluster")

test_that("cluster exercise works",{
  skip_on_cran()

  # import data
  data(ClusterExercise)

  ClusterExercise$Cluster.strat <- as.factor(ClusterExercise$Cluster.strat)
  # setup data
  # get rid of old region labels
  ClusterExercise$Area <- sum(unique(ClusterExercise$Area))
  ClusterExercise$Region.Label <- "Total"

  # fit detection function
  # df <- ds(ClusterExercise, truncation=1.5, key="hr",
  #          cutpoints=seq(0, 1.5, len=8), adjustment=NULL)

  # exactly as in Distance for Windows
  dat <- unflatten(ClusterExercise)
  dat <- dat$data[!is.na(dat$data$distance), ]
  dat <- dat[dat$distance<=1.5, ]
  dat <- Distance::create_bins(dat, seq(0, 1.5, len=8))
  # fit the exact function fitted by Distance
  result <- ddf(dsmodel = ~mcds(key = "hr", formula = ~1),
                data = dat, method = "ds",
                meta.data = list(width = 1.5, binned=TRUE,
                breaks=seq(0, 1.5, len=8)),
                control=list(initial=list(scale=log(0.7067555),
                                          shape=log(2.484188)),
                             nofit=TRUE))
  df <- result

  strat_N <- dht2(df, flatfile=ClusterExercise, strat_formula=~Cluster.strat,
                  stratification="object")

#            Estimate      %CV     df     95% Confidence Interval
#           ------------------------------------------------------
#  Stratum:  1
#  Hazard/Cosine
#    DS     0.11690E-01   26.38    36.13 0.69088E-02  0.19780E-01
#    D      0.11690E-01   26.38    36.13 0.69088E-02  0.19780E-01
#    N       8362.0       26.38    36.13  4942.0       14149.    
#  Stratum:  2
#  Hazard/Cosine
#    DS     0.98915E-02   32.15    31.43 0.52197E-02  0.18745E-01
#    D      0.19783E-01   32.15    31.43 0.10439E-01  0.37490E-01
#    N       14151.       32.15    31.43  7467.0       26817.    
#  Stratum:  3
#  Hazard/Cosine
#    DS     0.47959E-02   39.18    28.70 0.22132E-02  0.10392E-01
#    D      0.27876E-01   41.25    34.46 0.12462E-01  0.62358E-01
#    N       19940.       41.25    34.46  8914.0       44605.    
#  Pooled Estimates:
#            Estimate      %CV     df     95% Confidence Interval
#           ------------------------------------------------------
#    DS     0.26377E-01   20.41   117.85 0.17680E-01  0.39353E-01
#    D      0.59349E-01   24.51    76.50 0.36685E-01  0.96015E-01
#    N       42453.       24.51    76.50  26241.       68681.  
#
#            Estimate      %CV     df     95% Confidence Interval
#           ------------------------------------------------------
#  Stratum:  1
#    n       39.000
#    k       25.000
#    L       1842.8
#    n/L    0.21164E-01   23.72    24.00 0.13057E-01  0.34302E-01
#
#  Stratum:  2
#    n       33.000
#    k       25.000
#    L       1842.8
#    n/L    0.17908E-01   30.01    24.00 0.97691E-02  0.32826E-01
#
#  Stratum:  3
#    n       16.000
#    k       25.000
#    L       1842.8
#    n/L    0.86825E-02   37.45    24.00 0.41109E-02  0.18338E-01

  expect_equal(strat_N$Abundance, c(8362, 14151, 19940, 42453), tol=1e-2)
  expect_equal(strat_N$LCI, c(4942.0, 7467.0, 8914.0, 26241), tol=1e-2)
  expect_equal(strat_N$UCI, c(14149, 26817, 44605, 68681), tol=1e-2)
  expect_equal(strat_N$df, c(36.13, 31.43, 34.46, 76.50), tol=1e-1)
  expect_equal(strat_N$Abundance_CV, c(26.38, 32.15, 41.25, 24.51)/100,
               tol=1e-2)
})
