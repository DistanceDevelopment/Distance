# cue counting example

tol <- 1e-3
# make some results to check against
data(CueCountingExample)

# data fiddling
CueCountingExample$Effort <- CueCountingExample$Search.time

# format the cue data
cuedat <- CueCountingExample[,c("Cue.rate", "Cue.rate.SE", "Cue.rate.df")]
cuedat <- unique(cuedat)
names(cuedat) <- c("rate", "SE", "df")
# get rid of those columns, as we don't need them any more
CueCountingExample[, c("Cue.rate", "Cue.rate.SE", "Cue.rate.df",
                       "Sample.Fraction", "Sample.Fraction.SE")] <- list(NULL)
CueCountingExample$Label <- NULL

# set truncation
trunc <- 1.2

context("cue counting")

test_that("unbinned", {
  skip_on_cran()

#  Effort        :    36.47000    
#  # samples     :    92
#  Width         :    1.200000    
#  # observations:    40
# 
#  Model  1
#     Half-normal key, k(y) = Exp(-y**2/(2*A(1)**2))
# 
# 
#               Point        Standard    Percent Coef.        95% Percent
#   Parameter   Estimate       Error      of Variation     Confidence Interval
#   ---------  -----------  -----------  --------------  ----------------------
#     DS       0.81265E-01  0.31325E-01      38.55      0.36379E-01  0.18153    
#     E(S)      1.0000    
#     D        0.81265E-01  0.31325E-01      38.55      0.36379E-01  0.18153    
#     N         13653.       5262.8          38.55       6112.0       30498.    
#   ---------  -----------  -----------  --------------  ----------------------
# 
#  Measurement Units                
#  ---------------------------------
#  Density: Numbers/Sq. kilometers 
#      EDR: kilometers     
# 
#  Component Percentages of Var(D)
#  -------------------------------
#  Detection probability   :  21.7
#  Encounter rate          :  51.4
#  Cue rate                :  26.9


  # selected hn detection
  df_hn <- ds(CueCountingExample, truncation=trunc, key="hn", transect="point",
              adjustment=NULL, er_var="P3")

  # unstratified
  fs <- dht2(df_hn, flatfile=CueCountingExample,
             strat_formula=~1, multipliers=list(creation=cuedat),
             sample_fraction=0.5)
  expect_equal(fs$Abundance[nrow(fs)], 13653, tol=tol)
  expect_equal(fs$Abundance_se[nrow(fs)], 5262.8, tol=tol)
  expect_equal(fs$Abundance_CV[nrow(fs)], 38.55/100, tol=tol)
  expect_equal(fs$LCI[nrow(fs)], 6112.0, tol=tol)
  expect_equal(fs$UCI[nrow(fs)], 30498, tol=tol)



})


## compare with hr
#df_hr <- ds(CueCountingExample, truncation=trunc, key="hr", transect="point", adjustment=NULL)
#
#fs_st1_hr <- dht2(df_hr$ddf, flatfile=CueCountingExample,
#                  strat_formula=~1, multipliers=list(creation=cuedat),
#                  sample_fraction=0.5)
#fs_st1_hr


test_that("stratified estimation",{
  skip_on_cran()
#                         Estimate      %CV     df     95% Confidence Interval
#                        ------------------------------------------------------
# Stratum: 1. B                                              
# Half-normal/Cosine     
#                 DS     0.84002E-01   69.48    28.74 0.23261E-01  0.30335    
#                 D      0.84002E-01   69.48    28.74 0.23261E-01  0.30335    
#                 N       7140.0       69.48    28.74  1977.0       25785.    
# Stratum: 2. C                                              
# Half-normal/Cosine     
#                 DS     0.49123E-01   45.84    12.52 0.19051E-01  0.12667    
#                 D      0.49123E-01   45.84    12.52 0.19051E-01  0.12667    
#                 N       737.00       45.84    12.52  286.00       1900.0    
# Stratum: 3. D                                              
# Half-normal/Cosine     
#                 DS     0.94063E-01   44.06    15.37 0.38399E-01  0.23042    
#                 D      0.94063E-01   44.06    15.37 0.38399E-01  0.23042    
#                 N       1411.0       44.06    15.37  576.00       3456.0    
# Stratum: 4. E                                              
# Half-normal/Cosine     
#                 DS     0.84277E-01   60.38     7.10 0.22636E-01  0.31378    
#                 D      0.84277E-01   60.38     7.10 0.22636E-01  0.31378    
#                 N       1686.0       60.38     7.10  453.00       6276.0    
# Stratum: 5. F                                              
# Half-normal/Cosine     
#                 DS     0.78126E-01   63.80    29.43 0.23663E-01  0.25794    
#                 D      0.78126E-01   63.80    29.43 0.23663E-01  0.25794    
#                 N       2578.0       63.80    29.43  781.00       8512.0    
#
# Pooled Estimates:
#                         Estimate      %CV     df     95% Confidence Interval
#                        ------------------------------------------------------
#                 DS     0.80664E-01   45.23    50.57 0.33921E-01  0.19182    
#                 D      0.80664E-01   45.23    50.57 0.33921E-01  0.19182    
#                 N       13552.       45.23    50.57  5699.0       32226.    
#

  dat <- unflatten(CueCountingExample)
  dat <- dat$data
  # fit the exact function fitted by Distance
  result <- ddf(dsmodel = ~mcds(key = "hn", formula = ~1),
                data = dat, method = "ds",
                meta.data = list(width = 1.2, point=TRUE),
                control=list(initial=list(scale=log(0.4179)),
                             nofit=TRUE))
  df <- result

  fs_st1 <- dht2(df, flatfile=CueCountingExample,
                 strat_formula=~Region.Label, multipliers=list(creation=cuedat),
                 sample_fraction=0.5)

  expect_equal(fs_st1$Abundance,
               c(7140.0, 737.00, 1411.0, 1686.0, 2578.0, 13552), tol=tol)
  expect_equal(fs_st1$df, c(28.74, 12.52, 15.37, 7.10, 29.43, 50.57), tol=tol)
  expect_equal(fs_st1$Abundance_CV,
               c(69.48, 45.84, 44.06, 60.38, 63.80, 45.23)/100, tol=tol)
  expect_equal(fs_st1$LCI,
               c(1977.0, 286.00, 576.00, 453.00, 781.00, 5699.0), tol=tol)
  expect_equal(fs_st1$UCI,
               c(25785, 1900.0, 3456.0, 6276.0, 8512.0, 32226), tol=tol)

})
