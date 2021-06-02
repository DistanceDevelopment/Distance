# test ds()
par.tol <- 1e-5
N.tol <- 1e-3
lnl.tol <- 1e-4

context("summarize_ds_models()")

# data setup
data(book.tee.data)
egdata <- book.tee.data$book.tee.dataframe
egdata <- egdata[!duplicated(egdata$object),]

test_that("Error on different truncation distance",{
  set.seed(100)
  # some models
  t4 <- ds(egdata,4)
  t42 <- t4
  t4hr <- suppressWarnings(ds(egdata,4, key="hr"))
  t3 <- suppressWarnings(ds(egdata,3))
  t14 <- suppressWarnings(ds(egdata,list(left=1, right=4)))

  # when things are okay!
  out <- "           Model                                       Key function Formula\n1   \\\\texttt{t4}                                        Half-normal      ~1\n2  \\\\texttt{t42}                                        Half-normal      ~1\n3 \\\\texttt{t4hr} Hazard-rate with cosine adjustment term of order 2      ~1\n  C-vM p-value $\\\\hat{P_a}$ se($\\\\hat{P_a}$) $\\\\Delta$AIC\n1    0.7789695    0.5842744       0.04637627     0.000000\n2    0.7789695    0.5842744       0.04637627     0.000000\n3    0.8170780    0.5556728       0.07197512     2.501712"
  expect_output(print(summarize_ds_models(t4, t42, t4hr)), out, fixed=TRUE)

  # right truncation different
  expect_error(summarize_ds_models(t3,t4))
  # left
  expect_error(summarize_ds_models(t4,t14))
  # both
  expect_error(summarize_ds_models(t3,t4, t41))
})


test_that("Binning",{

  # largest cutpoint
  cp1 <- suppressMessages(ds(egdata,key="hn",order=0,
                  cutpoints=c(0,1,2,3,3.84)))
  cp11 <- cp1
  cp2 <- suppressMessages(ds(egdata,key="hn",order=0,
                  cutpoints=c(0,2,3,3.84)))

  out <- "           Model Key function Formula $\\\\chi^2$ $p$-value $\\\\hat{P_a}$\n1  \\\\texttt{cp1}  Half-normal      ~1           0.7452328    0.6340415\n2 \\\\texttt{cp11}  Half-normal      ~1           0.7452328    0.6340415\n  se($\\\\hat{P_a}$) $\\\\Delta$AIC\n1       0.05389734            0\n2       0.05389734            0"
  expect_output(print(summarize_ds_models(cp1, cp11)), out, fixed=TRUE)

  # different bins
  expect_error(summarize_ds_models(cp1,cp11,cp2))

  # mixing binned and unbinned
  ncp <- suppressMessages(ds(egdata,key="hn",order=0))
  expect_error(summarize_ds_models(cp1,ncp))

})
