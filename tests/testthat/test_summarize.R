# test ds()
par.tol <- 1e-5
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
  t4 <- ds(egdata, 4)
  t42 <- t4
  t4hr <- suppressWarnings(ds(egdata, 4, key="hr"))
  t3 <- suppressWarnings(ds(egdata, 3))
  t14 <- suppressWarnings(ds(egdata, list(left=1, right=4)))

  # when things are okay!
  out <- read.table(text="
\\texttt{t4}                                              Half-normal      ~1    0.7789695         0.5842744       0.04637627         0.000000
\\texttt{t42}                                             Half-normal      ~1    0.7789695         0.5842744       0.04637627         0.000000
\\texttt{t4hr} \"Hazard-rate with cosine adjustment term of order 2\"      ~1    0.8170780         0.5556728       0.07197555         2.501712
", header=FALSE)
  names(out) <- c("Model", "Key function", "Formula", "C-vM p-value",
                "$\\hat{P_a}$", "se($\\hat{P_a}$)", "$\\Delta$AIC")
  expect_equal(summarize_ds_models(t4, t42, t4hr), out, fixed=TRUE, tol=par.tol)

  # right truncation different
  expect_error(summarize_ds_models(t3, t4))
  # left
  expect_error(summarize_ds_models(t4, t14))
  # both
  expect_error(summarize_ds_models(t3, t4, t41))
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
  expect_equal(summarize_ds_models(cp1, cp11), out, fixed=TRUE, tol=par.tol)

  # different bins
  expect_error(summarize_ds_models(cp1, cp11, cp2))

  # mixing binned and unbinned
  ncp <- suppressMessages(ds(egdata, key="hn", order=0))
  expect_error(summarize_ds_models(cp1, ncp))

})
