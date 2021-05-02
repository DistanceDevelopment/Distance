# point transect exercise


# simulated data
data(PTExercise)
convert.units <- 0.01

context("PTexercise")

test_that("hn", {
# Effort        :    30.00000    
# # samples     :    30
# Width         :    20.00000    
# # observations:   131
#
# Model
#    Half-normal key, k(y) = Exp(-y**2/(2*A(1)**2))
#
#
#              Point        Standard    Percent Coef.        95% Percent
#  Parameter   Estimate       Error      of Variation     Confidence Interval
#  ---------  -----------  -----------  --------------  ----------------------
#    D         70.822       11.134          15.72       51.976       96.502    
#  ---------  -----------  -----------  --------------  ----------------------
#
# Measurement Units                
# ---------------------------------
# Density: Numbers/hectares       
#     EDR: meters         
#
# Component Percentages of Var(D)
# -------------------------------
# Detection probability   :  59.3
# Encounter rate          :  40.7
#
  # half-normal model
#PTExercise$size <- 1
#  df_hn <- ds(PTExercise, transect="point", key="hn", adjustment=NULL,
#              truncation=20, convert.units=convert.units)
#  df_hn <- dht2(df_hn, flatfile=PTExercise, strat_formula=~1)
#  expect_equal(attr(df_hn,"density")$Density, 70.822, tol=1e-1)
#  expect_equal(attr(df_hn,"density")$LCI, 51.976, tol=1e-2)
#  expect_equal(attr(df_hn,"density")$UCI, 96.502, tol=1e-2)
#  expect_equal(attr(df_hn,"density")$Density_se, 11.134, tol=1e-2)
#  expect_equal(attr(df_hn,"density")$Density_CV, .1572, tol=1e-2)

})


# hazard
test_that("hn", {
  df_hr <- ds(PTExercise, transect="point", key="hr", truncation=20,
              convert.units=convert.units, er.var="P3")
  # this gives answers per square metre
  df_hr$dht
})

test_that("hn", {
  df_unif <- ds(PTExercise, transect="point", key="unif", truncation=20,
                convert.units=convert.units, order=1, er.var="P3")
  # this gives answers per square metre
  df_unif$dht
})
