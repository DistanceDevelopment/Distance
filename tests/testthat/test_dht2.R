context("Testing dht2")

## first need to take the golf tee data
data(book.tee.data)
region <- book.tee.data$book.tee.region
samples <- book.tee.data$book.tee.samples
obs <- book.tee.data$book.tee.obs
obs <- obs[order(obs$object),]
egdata <- book.tee.data$book.tee.dataframe
# take only observer 1 data
egdata <- egdata[egdata$observer==1,]

obs <- merge(obs, egdata, by="object")
obs$Region.Label <- NULL
obs$distance <- NULL
obs$observer <- NULL
obs$detected <- NULL
obs$sex <- NULL
obs$exposure <- NULL

ds.dht.model <- suppressMessages(ds(egdata,4))

test_that("error thrown if geo stratum levels are mismatched",{

  region$Region.Label <- factor(region$Region.Label, 1:2, c("A","B"))
  expect_error(dht2(ds.dht.model, obs, samples, region,
                    strat_formula=~Region.Label),
    "Levels of Region.Label are different in `geo_strat` and `transects`")

})

test_that("undropped geo stratum covar levels don't cause errors",{

  region$sRegion.Label <- region$Region.Label
  samples$sRegion.Label <- samples$Region.Label

  region$Region.Label <- factor(region$sRegion.Label, 1:2, 1:2)
  samples$Region.Label <- factor(samples$sRegion.Label, 1:2, 1:2)
  d1 <- dht2(ds.dht.model, obs, samples, region, strat_formula=~Region.Label)

  region$Region.Label <- factor(region$sRegion.Label, 1:3, 1:3)
  samples$Region.Label <- factor(samples$sRegion.Label, 1:3, 1:3)

  d2 <- dht2(ds.dht.model, obs, samples, region, strat_formula=~Region.Label)

  expect_equivalent(d1, d2)

})


test_that("undropped individual stratum covar levels don't cause errors",{

  region <- data.frame(Area=sum(region$Area))

  # save sex var
  egdata$sexo <- egdata$sex

  # factorize
  egdata$sex <- factor(egdata$sexo, 0:1, 0:1)
  ds.dht.model <- suppressMessages(ds(egdata,4,formula=~sex))
  d1 <- dht2(ds.dht.model, obs, samples, region, strat_formula=~sex)

  egdata$sex <- factor(egdata$sexo, 0:2, 0:2)
  ds.dht.model <- suppressMessages(ds(egdata,4,formula=~sex))

  d2 <- dht2(ds.dht.model, obs, samples, region, strat_formula=~sex)

  expect_equivalent(d1, d2)
})

