# test ds()
par.tol <- 1e-4
N.tol <- 1e-3
lnl.tol <- 1e-4

context("summarize_ds_models()")

# data setup
data(book.tee.data)
egdata <- book.tee.data$book.tee.dataframe
egdata <- egdata[!duplicated(egdata$object), ]

test_that("Error on different truncation distance", {
  skip_on_cran()
  set.seed(100)
  # some models
  t4 <- ds(egdata, 4, optimizer = "R")
  t42 <- t4
  t4hr <- suppressWarnings(ds(egdata, 4, key="hr", optimizer = "R"))
  t3 <- suppressWarnings(ds(egdata, 3))
  t14 <- suppressWarnings(ds(egdata, list(left=1, right=4)))

  # when things are okay!
  out <- read.table(text="
\\texttt{t4}                                              Half-normal      ~1    0.7789695         0.5842744       0.04637627         0.000000
\\texttt{t42}                                             Half-normal      ~1    0.7789695         0.5842744       0.04637627         0.000000
\\texttt{t4hr} \"Hazard-rate with cosine adjustment term of order 2\"      ~1    0.8162846         0.5543305       0.07156778         2.454826
", header=FALSE)
  names(out) <- c("Model", "Key function", "Formula", "C-vM p-value",
                "$\\hat{P_a}$", "se($\\hat{P_a}$)", "$\\Delta$AIC")
  expect_equal(suppressWarnings(summarize_ds_models(t4, t42, t4hr)), out, fixed=TRUE, tol=par.tol)
  
  # right truncation different
  expect_error(suppressWarnings(summarize_ds_models(t3, t4)))
  # left
  expect_error(suppressWarnings(summarize_ds_models(t4, t14)))
  # both
  expect_error(suppressWarnings(summarize_ds_models(t3, t4, t41)))
})


test_that("Binning",{
  skip_on_cran()

  # largest cutpoint
  cp1 <- suppressMessages(ds(egdata, key="hn", order=0,
                  cutpoints=c(0, 1, 2, 3, 3.84)))
  cp11 <- cp1
  cp2 <- suppressMessages(ds(egdata, key="hn", order=0,
                  cutpoints=c(0, 2, 3, 3.84)))

  out <- read.table(text="
\\texttt{cp1}  Half-normal      ~1      0.8980817    0.6287014    0.05262633  0
\\texttt{cp11}  Half-normal     ~1      0.8980817    0.6287014    0.05262633  0
", header=FALSE)


  names(out) <- c("Model", "Key function", "Formula", "$\\chi^2$ $p$-value",
                  "$\\hat{P_a}$", "se($\\hat{P_a}$)", "$\\Delta$AIC")
  expect_equal(suppressWarnings(summarize_ds_models(cp1, cp11)), out, fixed=TRUE, tol=par.tol)

  # different bins
  expect_error(suppressWarnings(summarize_ds_models(cp1, cp11, cp2)))

  # mixing binned and unbinned
  ncp <- suppressMessages(ds(egdata, key="hn", order=0))
  expect_error(suppressWarnings(summarize_ds_models(cp1, ncp)))

})

test_that("Passing in models via a list",{
  skip_on_cran()
  data(book.tee.data)
  tee.data <- subset(book.tee.data$book.tee.dataframe, observer==1)
  ds.model <- ds(tee.data, 4)
  ds.model.cos <- ds(tee.data, 4, adjustment="cos", order=2)
  ds.model.hr <- ds(tee.data, 4, key = "hr", nadj = 0)
  
  expect_warning(test1 <- summarize_ds_models(ds.model, ds.model.cos, ds.model.hr), "Passing models via ... will be depricated in the next release, please pass models in a list using the models argument.")
  
  test2 <- summarize_ds_models(models = list(ds.model, ds.model.cos, ds.model.hr))
  
  expect_identical(test1[,2:7], test2[,2:7])
  
  test3 <- summarize_ds_models(list(ds.model, ds.model.cos, ds.model.hr))
  
  expect_identical(test1[,2:7], test3[,2:7])
  expect_identical(test3[,1], c("\\texttt{model 1}", "\\texttt{model 2}", "\\texttt{model 3}"))
  
  test4 <- summarize_ds_models(list(m1 = ds.model, m2 = ds.model.cos, m3 = ds.model.hr))
  expect_identical(test4[,1], c("\\texttt{m1}", "\\texttt{m2}", "\\texttt{m3}"))
})
