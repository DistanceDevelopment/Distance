context("Testing dht")

## first need to take the golf tee data
data(book.tee.data)
region <- book.tee.data$book.tee.region
samples <- book.tee.data$book.tee.samples
obs <- book.tee.data$book.tee.obs
obs <- obs[order(obs$object),]
egdata <- book.tee.data$book.tee.dataframe
# take only observer 1 data
egdata <- egdata[egdata$observer==1,]

test_that("dht without obs.table",{
  skip_on_cran()

  # first with obs table
  ds.dht.model<-suppressMessages(ds(egdata,4,region.table=region,
                                    monotonicity="strict",
                                    sample.table=samples,obs.table=obs))

  # check if you just forget obs.table an error is thrown
  # No longer throws error just displays to the user that it is only estimating
  # detection function
  #expect_error(suppressMessages(ds(egdata,4,region.table=region,
  #                                 monotonicity="strict",
  #                                 sample.table=samples)))

  # test that the result is the same if we do the correct merge
  dat <- merge(egdata, obs, by="object")
  ds.dht.model2 <- suppressMessages(ds(dat, 4, region.table=region,
                                       monotonicity="strict",
                                       sample.table=samples))
  expect_equivalent(ds.dht.model$dht, ds.dht.model2$dht)

})

# check density estimation when no area column present
# https://github.com/DistanceDevelopment/Distance/issues/56
test_that("no Area column works", {
  skip_on_cran()

  # with separate data.frames
  ds_with_Area <- ds(egdata, 4, region.table=region, adjustment=NULL,
                     sample.table=samples, obs.table=obs)

  region$Area <- NULL
  ds_no_Area <- ds(egdata, 4, region.table=region, adjustment=NULL,
                   sample.table=samples, obs.table=obs)

  # degrees of freedom/uncertainty will not match!
  # estimates should though!
  expect_equal(ds_with_Area$dht$individuals$D,#$Estimate,
               ds_no_Area$dht$individuals$D)#$Estimate)


  # with flatfile
  data(minke)

  minke_with_Area <- ds(minke, adjustment=NULL)

  minke$Area <- NULL
  minke_no_Area <- ds(minke, adjustment=NULL)

  # totals line doesn't match here as the proportion of effort between
  # strata is different
  expect_equal(minke_with_Area$dht$individuals$D[1:2,],
               minke_no_Area$dht$individuals$D[1:2,])

})
