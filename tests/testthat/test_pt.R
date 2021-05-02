# point transect testing

par.tol <- 1e-5
p.tol <- 1e-4
lnl.tol <- 1e-4

context("Point transects")

test_that("ptexample.distance from mrds gives same results",{

  set.seed(123)
  data(ptdata.distance)
  # mrds code for the same example
  #xx <- ddf(dsmodel=~cds(key="hn", formula=~1), data=ptdata.distance,
  #          method="ds",
  #          meta.data=list(point=TRUE,width=max(ptdata.distance$distance)))#,
  xx <- ds(ptdata.distance, key="hn", transect="point", adjustment=NULL)

  expect_that(xx$ddf$par,equals(2.283007,tol=par.tol))
  expect_that(xx$ddf$lnl,equals(-458.5701,tol=lnl.tol))
  expect_that(summary(xx$ddf)$average.p,equals(0.1644303,tol=p.tol))

})

test_that("ptexample.single from mrds gives same results -- binned",{

  set.seed(345)
  data(ptdata.single)

  # mrds code for the same example
  #ptdata.single$distbegin <- (as.numeric(cut(ptdata.single$distance,
  #                            10*(0:10)))-1)*10
  #ptdata.single$distend <- (as.numeric(cut(ptdata.single$distance,
  #                          10*(0:10))))*10
  #model <- ddf(data=ptdata.single, dsmodel=~cds(key="hn"),
  #                meta.data=list(point=TRUE,binned=TRUE,breaks=10*(0:10)))

  xx <-suppressMessages(ds(ptdata.single,key="hn",transect="point",
                           cutpoints=10*(0:10),order=0))

  expect_that(xx$ddf$par,equals(3.397266,tol=par.tol))
  expect_that(xx$ddf$lnl,equals(-687.7673,tol=lnl.tol))
  expect_that(summary(xx$ddf)$average.p,equals(0.1779359,tol=p.tol))

})

## amongst other things: test that left truncation does the right thing
## based on bug at https://github.com/DistanceDevelopment/Distance/issues/64
## this test is very slow so disabled for now
#test_that("camera trap data works",{
#
#  data(DuikerCameraTraps)
#
#  trunc.list <- list(left=2, right=15)
#  conversion <- convert_units("meter", NULL, "square kilometer")
#  DuikerCameraTraps.hr <- ds(DuikerCameraTraps, transect = "point",
#                             key="hr", adjustment = NULL,
#                             cutpoints = c(seq(2,8,1), 10, 12, 15),
#                             truncation = trunc.list, convert.units = conversion,
#                             er.var="P3")
#
#  mult.list <-list(creation=data.frame(rate=unique(DuikerCameraTraps$multiplier),
#                                       SE=0))
#  DuikerCameraTraps.hr.dens.mult <- dht2(DuikerCameraTraps.hr,
#                                         flatfile=DuikerCameraTraps,
#                                         strat_formula = ~1,
#                                         multipliers = mult.list,
#                                         convert_units = conversion)
#
#
#  #Check n, effort and number of samples
#  expect_equal(DuikerCameraTraps.hr.dens.mult$n, 5865)
#  expect_equal(DuikerCameraTraps.hr.dens.mult$Effort, 12317058)
#  expect_equal(DuikerCameraTraps.hr.dens.mult$k, 21)
#  # check density estimate
#  expect_equal(attr(DuikerCameraTraps.hr.dens.mult, "density")$Density,
#               14.51055, tol=lnl.tol)
#
#})
