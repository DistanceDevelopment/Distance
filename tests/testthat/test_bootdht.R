context("Testing bootdht")

test_that("area=0 with default summary",{

  skip_on_cran()

  data(amakihi)
  # subset for testing speed
  amakihi <- amakihi[1:300,]
  conv <- convert_units("meter", NULL, "hectare")
  surv <- ds(amakihi, transect="point", key="hr",
             formula=~Region.Label, convert.units=conv,
             truncation = 82.5, er.var="P3")
  expect_error(bootdht(surv, flatfile=amakihi, nboot=5),
               "No Area in flatfile, densities will be returned and the default summary function records only abundances. You need to write your own summary_fun.")

  amakihi$Area <- NULL
  expect_error(bootdht(surv, flatfile=amakihi, nboot=5),
               "No Area in flatfile, densities will be returned and the default summary function records only abundances. You need to write your own summary_fun.")

})


