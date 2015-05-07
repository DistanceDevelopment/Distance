context("Testing simplified data tables (tee data)")

## first need to take the golf tee data
data(book.tee.data)
region <- book.tee.data$book.tee.region
samples <- book.tee.data$book.tee.samples
obs <- book.tee.data$book.tee.obs
obs <- obs[order(obs$object),]
egdata <- book.tee.data$book.tee.dataframe
# take only observer 1 data
egdata <- egdata[egdata$observer==1,]

# merge sample onto region
flatdat <- merge(region,samples,by="Region.Label")
# merge obs table onto that
flatdat <- merge(flatdat,obs,by=c("Sample.Label","Region.Label"),all.x=TRUE)
# finally merge the distances onto that
flatdat <- merge(flatdat,egdata,by="object",all.x=TRUE)




test_that("Recover structure for tee data",{

  ## check we can recover the 4 data frames
  sepdat <- Distance:::checkdata(flatdat)

  # this should be a factor really, as is elsewhere
  obs$Region.Label <- as.factor(obs$Region.Label)

  ## now egdata should actually have the Region.Label and Sample.Label columns
  egdata <- merge(egdata, obs)
  egdata <- egdata[,c(1,9,8,2:7)]

  expect_equivalent(sepdat$region.table,region)
  expect_equivalent(sepdat$obs.table,obs)
  expect_equivalent(sepdat$sample.table,samples)
  expect_equivalent(sepdat$data,egdata)


})

test_that("Errors work",{

  flatdat[nrow(flatdat)+1,] <- flatdat[1,]
  flatdat[nrow(flatdat),]$Effort <- 1322

  expect_error(Distance:::checkdata(flatdat))

  # check we get an error if the columns are right but their case is wrong
  flatdat$effort <- flatdat$Effort
  flatdat$Effort <- NULL
  expect_error(Distance:::checkdata(flatdat))

})
#test_that("Extra regions",{
#
#  # add an extra region
#  region.e <- data.frame(Region.Label=1:3,Area=c(1040,640,500))
#  samples
#
#  # estimate using all tables
#  dsmod.e <- ds(egdata,4,region.table=region.e, sample.table=samples,
#                obs.table=obs)
#
#}
