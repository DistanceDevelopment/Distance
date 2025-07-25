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

# make a flatfile version
obs2 <- obs
obs2$size <- NULL
ff <- merge(region, samples, by="Region.Label", all.x=TRUE)
ff <- merge(ff, obs2, by="Sample.Label", all.x=TRUE)
ff <- merge(ff, egdata, by="object", all.x=TRUE)

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
#ff <- ff[, c("Region.Label", "Area","Sample.Label", "Effort", "size",
#             "distance")]
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

  ff_formfac <- dht2(ds.dht.model, flatfile=ff,
                  strat_formula=~as.factor(Region.Label))
  formfac <- dht2(ds.dht.model, obs, samples, region,
                  strat_formula=~as.factor(Region.Label))


  ff$Region.Label <- as.factor(ff$Region.Label)
  ff_prefac <- dht2(ds.dht.model, flatfile=ff, strat_formula=~Region.Label)
  prefac <- dht2(ds.dht.model, obs, samples, region,
                 strat_formula=~Region.Label)

  expect_equivalent(prefac, formfac)
  expect_equivalent(ff_prefac, ff_formfac)
})

test_that("error thrown if stratum name doesn't exist in flatfile",{

  expect_error(dht2(ds.dht.model, flatfile=ff,
                    strat_formula=~boop),
               "Column\\(s\\): boop not in \\`flatfile\\`")

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


# use amakihi for the next few examples
data(amakihi)
conv.am <- convert_units("meter", NULL, "hectare")
easy.am <- ds(amakihi, transect="point", key="hr", convert_units = conv.am,
              adjustment=NULL, er_var="P3")
uf <- unflatten(amakihi)

test_that("density estimation, no innes", {

  # dht
  dht_old <- dht(easy.am$ddf, uf$region.table, uf$sample.table,
                 options=list(convert.units=conv.am, ervar="P3", varflag=1))
  dht2_zero <- dht2(easy.am, flatfile=amakihi, convert_units = conv.am,
                    strat_formula = ~Region.Label, innes=FALSE, er_est="P3")

  expect_equal(dht_old$individuals$D$Estimate, dht2_zero$Abundance)
  expect_equal(dht_old$individuals$summary$ER, dht2_zero$ER)
  expect_equal(dht_old$individuals$summary$n, dht2_zero$n)
  expect_equal(dht_old$individuals$summary$CoveredArea, dht2_zero$Covered_area)
  #expect_equal(dht_old$individuals$summary$se.ER, sqrt(dht2_zero$ER_var))
  expect_equal(dht_old$individuals$D$se, dht2_zero$Abundance_se, tol=1e-4)

})

test_that("density estimation, no innes, longform", {

  # dht
  dht_old <- dht(easy.am$ddf, uf$region.table, uf$sample.table,
                 options=list(convert.units=conv.am, ervar="P3", varflag=1))
  # now with dht2 equiv
  dht2_zero <- dht2(easy.am, observations=uf$obs.table,
                    transects=uf$sample.table, uf$region.table,
                    convert_units = conv.am,
                    strat_formula = ~Region.Label, innes=FALSE, er_est="P2")

  expect_equal(dht_old$individuals$D$Estimate, dht2_zero$Abundance)
  expect_equal(dht_old$individuals$summary$ER, dht2_zero$ER)
  expect_equal(dht_old$individuals$summary$n, dht2_zero$n)
  expect_equal(dht_old$individuals$summary$CoveredArea, dht2_zero$Covered_area)
  #expect_equal(dht_old$individuals$summary$se.ER, sqrt(dht2_zero$ER_var))
#  expect_equal(dht_old$individuals$summary$cv.ER, sqrt(dht2_zero$ER_CV))
  expect_equal(dht_old$individuals$D$se, dht2_zero$Abundance_se, tol=1e-3)

})


test_that("density estimation, innes", {

  # dht
  dht_old <- dht(easy.am$ddf, uf$region.table, uf$sample.table, uf$obs.table,
                 options=list(convert.units=conv.am, ervar="P2", varflag=2))
  dht2_zero <- dht2(easy.am, flatfile=amakihi, convert_units = conv.am,
                    strat_formula = ~Region.Label, innes=TRUE, er_est="P2")

  expect_equal(dht_old$individuals$D$Estimate, dht2_zero$Abundance)
  expect_equal(dht_old$individuals$summary$ER, dht2_zero$ER)
  expect_equal(dht_old$individuals$summary$CoveredArea, dht2_zero$Covered_area)
  #expect_equal(dht_old$individuals$summary$se.ER, sqrt(dht2_zero$ER_var))
  expect_equal(dht_old$individuals$D$se, dht2_zero$Abundance_se, tol=1e-5)

})

test_that("varflag=0 works", {
  # dht with varflag=0
  dht_old <- dht(easy.am$ddf, uf$region.table, uf$sample.table, uf$obs.table,
                 options=list(varflag=0, convert.units=conv.am, ervar="P3"))
  # now with dht2 equiv
  dht2_zero <- dht2(easy.am, flatfile=amakihi, convert_units = conv.am,
                    strat_formula = ~Region.Label, binomial_var=TRUE)

  expect_equal(dht_old$individuals$D$Estimate, dht2_zero$Abundance)
  expect_equal(dht_old$individuals$D$se, dht2_zero$Abundance_se, tol=1e-5)
  expect_equal(dht_old$individuals$D$cv, dht2_zero$Abundance_CV, tol=1e-5)
  expect_equal(dht_old$individuals$D$df, dht2_zero$df)
  expect_equal(dht_old$individuals$D$lcl, dht2_zero$LCI, tol=1e-5)
  expect_equal(dht_old$individuals$D$ucl, dht2_zero$UCI, tol=1e-5)
})


data(minke)

test_that("0 effort errors", {

  minke$object <- NA
  minke$object[!is.na(minke$distance)] <- 1:sum(!is.na(minke$distance))

  minke_df <- ds(minke, truncation=1.5, adjustment=NULL)
  minke_dht2 <- dht2(minke_df, flatfile=minke, stratification="geographical",
                     strat_formula=~Region.Label)
  minke$Effort[minke$Region.Label=="South"] <- 0

  # if some effort is actually zero
  expect_error(dht2(minke_df, flatfile=minke, stratification="geographical",
                    strat_formula=~Region.Label),
               "Some transects have Effort <=0 or NA")

})

test_that("can accept data.frame and scalar sampling fractions", {
  
  # Add a sampling fraction
  minke.samp.frac <- data.frame(Sample.Label = unique(minke$Sample.Label))
  minke.samp.frac$fraction <- c(rep(0.5,15),rep(1,10))
  # first fitting the detection function
  minke_df <- ds(minke, truncation=1.5, adjustment=NULL)
  # now estimate abundance using dht2
  # stratum labels are in the Region.Label column
  minke_dht2 <- dht2(minke_df, 
                     flatfile=minke, 
                     stratification="geographical", 
                     sample_fraction = minke.samp.frac,
                     strat_formula=~Region.Label)
  # could compare this to minke_df$dht and see the same results
  minke_dht2
  # Alternative
  minke_dht2 <- dht2(minke_df, 
                     flatfile=minke, 
                     stratification="geographical", 
                     sample_fraction = 0.5,
                     strat_formula=~Region.Label)
  # could compare this to minke_df$dht and see the same results
  minke_dht2
  
})


test_that("Deal with strata with no detections", {
  
  # Add a third Region with two samples but no observations.
  # Adding two samples instead of just one to prevent a warning unstable ER var
  new.row1 <- minke[1, ] %>%
    mutate(
      Region.Label = "East",
      Area = 10000,
      Sample.Label = max(minke$Sample.Label) + 1,
      Effort = 10,
      distance = NA,
      object = NA
    )
  
  new.row2 <- minke[1, ] %>%
    mutate(
      Region.Label = "East",
      Area = 10000,
      Sample.Label = max(minke$Sample.Label) + 2,
      Effort = 10,
      distance = NA,
      object = NA
    )
  
  minke <- rbind(minke, new.row1, new.row2)
  # Fit a detection function
  detfn <- ds(data=minke, truncation=1.5, key="hr", adjustment=NULL)
  
  result <- dht2(ddf = detfn,
                 flatfile = minke,
                 strat_formula =  ~ Region.Label,
                 stratification = "geographical")
  
  expect_true(is(result, "dht_result"))
  
  # Check values in output
  new.row3 <- minke[1, ] %>%
    mutate(
      Region.Label = "East",
      Area = 10000,
      Sample.Label = max(minke$Sample.Label) + 3,
      Effort = 10,
      distance = NA,
      object = NA
    )
  
  minke2 <- rbind(minke, new.row1, new.row2, new.row3)
  detfn <- ds(data=minke2, truncation=1.5, key="hr", adjustment=NULL, optimizer = "R")
  dht2.results <- dht2(ddf = detfn,
                       flatfile = minke2,
                       strat_formula =  ~ Region.Label,
                       stratification = "geographical")
  
  expect_equal(dht2.results$Abundance[1], 0)
  expect_equal(dht2.results$Abundance_se[4], 4834.4118, , tolerance = 0.0001)
  expect_equal(dht2.results$df[4],15.25496, tolerance = 0.00001)
  
})
