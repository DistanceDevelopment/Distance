# cue counting example

# test numerical tollerance
tol <- 1e-6

# make some results to check against
data(golftees)

# data fiddling
gtees <- golftees[golftees$observer==1 & golftees$detected==1, ]
#gtees$size <- 1
dat <- unflatten(gtees)

trunc <- 4

make_old_summ_cluster <- function(object){

  # print out as dht
  object$Region <- as.factor(object$Region)
  object$CoveredArea <- object$Covered_area
  object$se.ER <- sqrt(object$ER_var)
  object$cv.ER <- object$ER_CV
  object$se.mean <- sqrt(object$group_var)
  object$mean.size <- object$group_mean

  class(object) <- "data.frame"

  summ <- object[, c("Region", "Area", "CoveredArea", "Effort", "n",
                     "ER", "se.ER", "cv.ER", "mean.size", "se.mean")]
  summ$se.ER[is.na(summ$se.ER) | is.nan(summ$se.ER)] <- 0
  summ$cv.ER[is.na(summ$cv.ER) | is.nan(summ$cv.ER)] <- 0
  return(summ)
}
make_old_abund_individual <- function(object){

  object$Label <- as.factor(object$Region)
  object$Region <- NULL
  object$Estimate <- object$Abundance
  object$cv <- object$Abundance_CV
  object$se <- object$Abundance_se
  object$lcl <- object$LCI
  object$ucl <- object$UCI
  object$df <- object$df
  class(object) <- "data.frame"
  object[, c("Label", "Estimate", "se", "cv", "lcl", "ucl", "df")]
}


context("golftees")

test_that("ER variance", {

  # output using old dht
  df <- ds(gtees, truncation=trunc, key="hn", adjustment=NULL)

  # now do a fancy thing
  dat$obs.table <- dat$obs.table[dat$obs.table$object %in% gtees$object, ]
  fs_st1 <- expect_warning(dht2(df$ddf, dat$obs.table, dat$sample.table,
                                dat$region.table, strat_formula=~Region.Label,
                                innes=FALSE),
                           "One or more strata have only one transect, cannot calculate empirical encounter rate variance")

  # tests
  oldres <- df$dht$clusters$summary
  oldres$k <- NULL
  oldres$mean.size <- df$dht$individuals$summary$mean.size
  oldres$se.mean <- df$dht$individuals$summary$se.mean
  expect_equal(oldres,
               make_old_summ_cluster(fs_st1), tolerance=tol)

  fs_st1$Region <- "Total"
  # TODO: this is very stupid and probably a bug in mrds
#  aa <- df$dht$individuals$N
#  aa[, -1] <- as.numeric(aa[, -1])
#  expect_equal(aa, make_old_abund_individual(fs_st1),
#               tolerance=tol, check.attributes=FALSE)

})


test_that("Same results as Distance", {
  gtees$sex <- as.factor(gtees$sex)
  gtees$sex <- relevel(gtees$sex, ref="1")
  df <- ds(gtees, truncation=trunc, key="hn", adjustment=NULL, formula=~sex)

  fs_st1 <- expect_warning(dht2(df$ddf, dat$obs.table, dat$sample.table,
                                dat$region.table, strat_formula=~Region.Label,
                                innes=FALSE),
                           "One or more strata have only one transect, cannot calculate empirical encounter rate variance")



  df <- ds(gtees, truncation=trunc, key="hn", adjustment=NULL, formula=~sex,
           initial_values=list(scale=c(log(2.30706), -0.583378)))

    result <- ddf(dsmodel = ~mcds(key = "hn", formula = ~sex, adj.scale="width"),
                   data = gtees, method = "ds",
                   meta.data = list(width = trunc),
                   control=list(initial=list(scale=c(log(2.306594), -0.5832343)), nofit=TRUE))





# redo DISTANCE's calculations
devtools::load_all("~/current/mrds")
par <- result$par
fl <- function(x, ...){
  ddfobj <- assign.par(list(...)$ddfobj, x)
  -flpt.lnl(x, ddfobj=ddfobj, misc.options=list(...)$misc.options)
}

n <- nrow(result$ds$aux$ddfobj$xmat)

misc<-list(point=result$meta.data$point, int.range=result$meta.data$int.range,
           showit=result$control$showit,
           integral.numeric=result$control$integral.numeric,
           breaks=result$meta.data$breaks, maxiter=result$control$maxiter,
           refit=result$control$refit, nrefits=result$control$nrefits,
           parscale=result$control$parscale, mono=result$meta.data$mono,
           mono.strict=result$meta.data$mono.strict, binned=result$meta.data$binned,
           width=result$meta.data$width, standardize=result$control$standardize,
           mono.points=result$control$mono.points, mono.tol=result$control$mono.tol,
           mono.delta=result$control$mono.delta, debug=result$control$debug,
           nofit=result$control$nofit, left=result$meta.data$left,
           silent=result$control$silent,
           mono.random.start=result$control$mono.random.start
          )



jacob <- numDeriv::jacobian(fl, x=par, ddfobj=result$ds$aux$ddfobj, misc.options=misc)

H <- (t(jacob) %*% jacob)/n

Hinv <- mrds:::solvecov(H)$inv


#se <- sqrt(diag(Hinv)/n)
#se



# make a sandwich
f <- function(x) c(exp(x[1]), x[2])

bread <- numDeriv::jacobian(f, x=par)


t(bread) %*% Hinv %*% bread


sqrt(diag(t(bread) %*% Hinv %*% bread)/n)

z <- result$ds$aux$ddfobj$scale$dm

dfun <- function(x, z, par){
  sigma <- z%*%par
  D <- (x/sigma^2)/2
  valkey <- exp(D)
  DKEYDS <- x^2 * valkey/sigma^2

  return(DKEYDS)
}

dfun(result$data$distance, z, result$par)



# corr to cov

mm <- as.matrix(read.table(text="
1 -0.75969
-0.75969 1
"))

D <- diag(se)

# distance vcov
D%*%mm%*%D

# send to Len
#f0 <- function(x, z, ddfobj){
#  key.scale <- exp(z%*%x)
#  g <- keyfct.hn(0, key.scale)
#  ddfobj <- assign.par(ddfobj, x)
#  mu <- integratepdf(ddfobj, NULL, trunc, c(0, trunc), standardize=FALSE)
#  g/mu
#}
f0 <- function(x, model){
  model$par <- x
  model$ds$aux$ddfobj <- assign.par(model$ds$aux$ddfobj, x)
  1/predict(model, esw=TRUE, compute=TRUE)$fitted
}

f0(result$par, result)


bread_f0 <- (numDeriv::jacobian(f0, x=result$par, model=result))[c(1,3), ]


t(bread_f0) %*% (Hinv/n) %*% bread_f0


bread_f0 <- (numDeriv::jacobian(f0, x=result$par, model=result))[1, ]
bread_f0_2 <- (numDeriv::jacobian(f0, x=result$par, model=result))[3, ]


vf0 <- sum(c(
rep((t(bread_f0) %*% (Hinv/n) %*% bread_f0)[1,1], sum(gtees$sex==1)),
rep((t(bread_f0_2) %*% (Hinv/n) %*% bread_f0_2)[1,1], sum(gtees$sex==0))))/124^2


f0a <- function(x, n, model){
  model$par <- x
  model$ds$aux$ddfobj <- assign.par(model$ds$aux$ddfobj, x)
  Nhat <- sum(1/predict(model, compute=TRUE)$fitted)
  1/(n/Nhat * 4)
}


bread_f0 <- (numDeriv::jacobian(f0a, x=result$par, model=result, n=124))

(bread_f0 %*% (Hinv/n) %*% t(bread_f0))

sqrt(bread_f0 %*% (Hinv/n) %*% t(bread_f0))/f0a(result$par, 124, result)



pa <- function(x, n, model){
  model$par <- x
  model$ds$aux$ddfobj <- assign.par(model$ds$aux$ddfobj, x)
  Nhat <- sum(1/predict(model, compute=TRUE)$fitted)
  n/Nhat
}


bread_p <- (numDeriv::jacobian(pa, x=result$par, model=result, n=n))

(bread_p %*% (Hinv/n) %*% t(bread_p))

sqrt(bread_p %*% (Hinv/n) %*% t(bread_p))/pa(result$par, n, result)




p <- function(x, n, model){
  model$par <- x
  model$ds$aux$ddfobj <- assign.par(model$ds$aux$ddfobj, x)
  predict(model, compute=TRUE)$fitted
}


bread_p <- (numDeriv::jacobian(p, x=result$par, model=result))[1, ]
bread_p_2 <- (numDeriv::jacobian(p, x=result$par, model=result))[3, ]


vp <- sum(c(
rep((t(bread_p) %*% (Hinv/n) %*% bread_p)[1,1], sum(gtees$sex==1)),
rep((t(bread_p_2) %*% (Hinv/n) %*% bread_p_2)[1,1], sum(gtees$sex==0))))/124^2

sqrt(vp)/pa(result$par, n, result)



