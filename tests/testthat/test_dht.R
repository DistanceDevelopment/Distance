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

  # first with obs table
  ds.dht.model<-suppressMessages(ds(egdata,4,region.table=region, monotonicity="strict", sample.table=samples,obs.table=obs))

  # check if you just forget obs.table an error is thrown
  # No longer throws error just displays to the user that it is only estimating detection function
  #expect_error(suppressMessages(ds(egdata,4,region.table=region, monotonicity="strict", sample.table=samples)))

  # test that the result is the same if we do the correct merge
  dat <- merge(egdata, obs, by="object")
  ds.dht.model2 <- suppressMessages(ds(dat, 4, region.table=region,
                                       monotonicity="strict",
                                       sample.table=samples))
  expect_equivalent(ds.dht.model$dht, ds.dht.model2$dht)

})

