# test if sample fractions work

context("Sample fractions")

data(wren_lt)

ll <- ds(wren_lt, convert_units=convert_units("meter", "kilometer", "hectare"),
         adjustment=NULL)

test_that("single sample fractions",{

  samfrac <- dht2(ll, flatfile=wren_lt, strat_formula=~1, sample_fraction=0.5,
                  convert_units=convert_units("meter", "kilometer", "hectare"))

  wren_lt$Effort <- wren_lt$Effort * 0.5
  manual <- dht2(ll, flatfile=wren_lt, strat_formula=~1,
                  convert_units=convert_units("meter", "kilometer", "hectare"))

  expect_equal(manual$Covered_area, samfrac$Covered_area)
  expect_equal(manual$Abundance_se, samfrac$Abundance_se)
  expect_equal(manual$Abundance, samfrac$Abundance)

})


test_that("data.frame sample fractions",{

  sf_df <- data.frame(Sample.Label = unique(wren_lt$Sample.Label),
                      fraction     =  0.5)

  samfrac <- dht2(ll, flatfile=wren_lt, strat_formula=~1, sample_fraction=sf_df,
                  convert_units=convert_units("meter", "kilometer", "hectare"))

  wren_lt$Effort <- wren_lt$Effort * 0.5
  manual <- dht2(ll, flatfile=wren_lt, strat_formula=~1,
                  convert_units=convert_units("meter", "kilometer", "hectare"))

  expect_equal(manual$Covered_area, samfrac$Covered_area)
  expect_equal(manual$Abundance_se, samfrac$Abundance_se)
  expect_equal(manual$Abundance, samfrac$Abundance)

})

test_that("errors",{

  sf_df <- data.frame(Sample.Label = unique(wren_lt$Sample.Label),
                      friction     =  0.5)

  expect_error(samfrac <- dht2(ll, flatfile=wren_lt, strat_formula=~1,
                               sample_fraction=sf_df,
                               convert_units=convert_units("meter", "kilometer",
                                                           "hectare")),
               "sample_fraction data.frame columns must be \"Sample.Label\" and \"fraction\"")


  expect_error(samfrac <- dht2(ll, flatfile=wren_lt, strat_formula=~1,
                               sample_fraction=c(0.1, 0.4),
                               convert_units=convert_units("meter", "kilometer",
                                                           "hectare")),
               "sample_fraction must be a single number or a data.frame")

  expect_error(samfrac <- dht2(ll, flatfile=wren_lt, strat_formula=~1,
                               sample_fraction=-0.1,
                               convert_units=convert_units("meter", "kilometer",
                                                           "hectare")),
               "sample_fraction must be > 0")

})
