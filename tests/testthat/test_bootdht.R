context("Testing bootdht")

test_that("area=0 with default summary",{

  skip_on_cran()

  data(amakihi)
  # subset for testing speed
  amakihi <- amakihi[1:300,]
  conv <- convert_units("meter", NULL, "hectare")
  surv <- ds(amakihi, transect="point", key="hr",
             formula=~Region.Label, convert_units=conv,
             truncation = 82.5, er_var="P3")
  expect_error(bootdht(surv, flatfile=amakihi, nboot=5),
               "No Area in flatfile, densities will be returned and the default summary function records only abundances. You need to write your own summary_fun.")

  amakihi$Area <- NULL
  expect_error(bootdht(surv, flatfile=amakihi, nboot=5),
               "No Area in flatfile, densities will be returned and the default summary function records only abundances. You need to write your own summary_fun.")

})


# generate some data to test on
dat <- data.frame(Sample.Label = 1:10)
dat$Region.Label <- c(rep(1, 5), rep(2, 5))
set.seed(123)
dat2 <- dat[sample(1:nrow(dat), nrow(dat), replace=TRUE), ]
dat2$object <- 1:nrow(dat2)
dat <- dat[!(dat$Sample.Label %in% dat2$Sample.Label),]
dat$object <- NA
dat <- rbind(dat, dat2)
dat$save.ind <- 1:nrow(dat)
dat$distance <- rnorm(nrow(dat), 0, 50)
dat$distance <- ifelse(is.na(dat$object),NA,dat$distance)


test_that("resamples work - strata", {

  set.seed(123)
  obs <- bootdht_resample_data(dat, c("Region.Label"))

  expect_equal(obs$save.ind, c(1, 4, 5, 7, 9, 10, 2, 3, 6, 13, 8, 11, 12))
})


test_that("resamples work - sample", {

  set.seed(123)
  obs <- bootdht_resample_data(dat, c("Sample.Label"))

  expect_equal(obs$save.ind, c(7, 7, 4, 5, 4, 5, 7, 12, 8, 11, 2, 3, 6, 13))
})

test_that("resamples work - object", {

  set.seed(123)
  obs <- bootdht_resample_data(dat, c("object"))

  expect_equal(obs$save.ind, c(1, 4, 5, 6, 5, 2, 3, 6, 11, 9))
})

test_that("resamples work - strata/object", {

  set.seed(123)
  obs <- bootdht_resample_data(dat, c("Region.Label", "object"))

  expect_equal(obs$save.ind, c(1, 4, 5, 5, 5, 9, 2, 3, 13, 8, 11, 8))
})


test_that("resamples work - sample/object", {

  set.seed(123)
  obs <- bootdht_resample_data(dat, c("Sample.Label", "object"))

  expect_equal(obs$save.ind, c(7, 4, 4, 4, 5, 4, 12, 8, 2, 3, 6, 13))
})

test_that("resamples work - strata/sample/object", {

  set.seed(1223)
  obs <- bootdht_resample_data(dat, c("Region.Label", "Sample.Label", "object"))

  expect_equal(obs$save.ind, c(6, 1, 9, 6, 6, 2, 3, 8, 11, 6, 13, 8))
})

# generate some data to test on
dat <- data.frame(Sample.Label = 1:10)
dat$Region.Label <- c(rep(1, 5), rep(2, 5))
set.seed(115)
dat2 <- dat
dat2$Region.Label <- paste("YEAR", as.character(dat2$Region.Label), sep = "")
dat2 <- rbind(dat2, dat2, dat2)
dat2$distance <- abs(rnorm(nrow(dat2), 0, 25))
dat2$Effort <- 1000
dat2$Area <- 100
dat2$size <- rpois(nrow(dat2), 20)
dat2$object <- 1:nrow(dat2)

conversion.factor <- convert_units("meter", "meter", "square kilometer")

fit.ds <- ds(data=dat2,
             truncation=50,
             key="hn",
             adjustment=NULL,
             convert_units=conversion.factor,
             formula=~size)



test_that("Issue #158 is fixed (stratum names > 'Total' bug)", {
  
  skip_on_cran()
  
  set.seed(225)
  easybootn <- suppressMessages(bootdht(model=fit.ds,
                       flatfile=dat2,
                       summary_fun=bootdht_Nhat_summarize,
                       convert_units=conversion.factor,
                       sample_fraction=1,
                       nboot=3, cores=1,
                       progress_bar = "none"))

  # Make check table
  check.tab <- easybootn %>% dplyr::group_by(Label) %>%
    dplyr::summarize(LCI=quantile(Nhat,probs=0.025,na.rm=TRUE),
                     UCI=quantile(Nhat,probs=0.975,na.rm=TRUE))
  
  # This was TRUE should now be FALSE after fix
  expect_false(any(is.na(check.tab$LCI)))
  
  # Now check cluster fix
  bootdht_NChat_summarize <- function(ests, fit) {
    return(data.frame(Label = ests$clusters$N$Label,
                      Nhat  = ests$cluster$N$Estimate))
  }
  
  set.seed(225)
  easybootn <- suppressMessages(bootdht(model=fit.ds,
                       flatfile=dat2,
                       summary_fun=bootdht_NChat_summarize,
                       convert_units=conversion.factor,
                       sample_fraction=1,
                       nboot=3, cores=1,
                       progress_bar = "none"))
  # Make check table
  check.tab <- easybootn %>% dplyr::group_by(Label) %>%
    dplyr::summarize(LCI=quantile(Nhat,probs=0.025,na.rm=TRUE),
                     UCI=quantile(Nhat,probs=0.975,na.rm=TRUE))
  
  # This was TRUE should now be FALSE after fix
  expect_false(any(is.na(check.tab$LCI)))
})


# Test data from Eric
data(minke)
# convert exact distances into bins
vals <- seq(0,2,.2)
minke$bin <- cut(minke$distance, breaks=seq(0, 2, .2), right=FALSE, labels=FALSE)
minke$distbegin <- vals[minke$bin]
minke$distend <- vals[minke$bin+1]
# remove exact distances
minke$distance <- NULL

test_that("Issue #158 is fixed (stratum names > 'Total' bug)", {
  
  skip_on_cran()
  mod1 <- ds(minke)
  set.seed(225)
  bootout <- bootdht(mod1, flatfile=minke,  nboot=3)
  expect_true(nrow(bootout) > 0)
})