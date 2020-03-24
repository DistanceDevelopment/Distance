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

obsold <- obs
obs <- merge(obs, egdata, by="object")
obs$Region.Label <- NULL
obs$distance <- NULL
obs$observer <- NULL
obs$detected <- NULL
obs$sex <- NULL
obs$exposure <- NULL

ds.dht.model <- suppressMessages(ds(egdata,4))


#test_that("error thrown if geo strat but object",{
#
#
#obs2 <- obs
#obs2$size <- NULL
#ff <- merge(region, samples, by="Region.Label", all.x=TRUE)
#ff <- merge(ff, obs2, by="Sample.Label", all.x=TRUE)
#ff <- merge(ff, egdata, by="object", all.x=TRUE)
#
#ff <- ff[, c("Region.Label", "Area","Sample.Label", "Effort", "size", "distance")]
#
#dd <- dht2(ds.dht.model, flatfile=ff, strat_formula=~Region.Label)
#
#dd <- dht2(ds.dht.model, flatfile=ff,
#           strat_formula=~Region.Label, stratification="object")
#
#dd <- dht2(ds.dht.model, flatfile=ff,
#           strat_formula=~as.factor(Region.Label))
#
#ff$fa <- as.factor(ff$Region.Label)
#dd <- dht2(ds.dht.model, flatfile=ff,
#           strat_formula=~fa)
#
##  expect_error(dht2(ds.dht.model, obs, samples, region,
##                    strat_formula=~Region.Label, stratification="object"),
##    "Levels of Region.Label are different in `geo_strat` and `transects`")
#
#  dd <- dht2(ds.dht.model, obs, samples, region,
#             strat_formula=~as.factor(Region.Label))
#})

test_that("as.factor works in formula",{

  obs2 <- obs
  obs2$size <- NULL
  ff <- merge(region, samples, by="Region.Label", all.x=TRUE)
  ff <- merge(ff, obs2, by="Sample.Label", all.x=TRUE)
  ff <- merge(ff, egdata, by="object", all.x=TRUE)
  ff_formfac <- dht2(ds.dht.model, flatfile=ff,
                  strat_formula=~as.factor(Region.Label))
  formfac <- dht2(ds.dht.model, obs, samples, region,
                  strat_formula=~as.factor(Region.Label))


  ff$Region.Label <- as.factor(ff$Region.Label)
  ff_prefac <- dht2(ds.dht.model, flatfile=ff, strat_formula=~Region.Label)
  prefac <- dht2(ds.dht.model, obs, samples, region, strat_formula=~Region.Label)

  expect_equivalent(prefac, formfac)
  expect_equivalent(ff_prefac, ff_formfac)
})


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





