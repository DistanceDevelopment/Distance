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

  data(ptdata.single)

  set.seed(345)

  # mrds code for the same example
  #ptdata.single$distbegin <- (as.numeric(cut(ptdata.single$distance,10*(0:10)))-1)*10
  #ptdata.single$distend <- (as.numeric(cut(ptdata.single$distance,10*(0:10))))*10
  #model <- ddf(data=ptdata.single, dsmodel=~cds(key="hn"),
  #                meta.data=list(point=TRUE,binned=TRUE,breaks=10*(0:10)))

  xx <-suppressMessages(ds(ptdata.single,key="hn",transect="point",
                           cutpoints=10*(0:10),order=0))

  expect_that(xx$ddf$par,equals(3.397266,tol=par.tol))
  expect_that(xx$ddf$lnl,equals(-687.7673,tol=lnl.tol))
  expect_that(summary(xx$ddf)$average.p,equals(0.1779359,tol=p.tol))

})
