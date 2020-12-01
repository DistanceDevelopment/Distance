`Distance`
==========

[![Build Status](https://travis-ci.org/DistanceDevelopment/Distance.svg?branch=master)](https://travis-ci.org/DistanceDevelopment/Distance)
[![CRAN (RStudio Mirror) Downloads](http://cranlogs.r-pkg.org/badges/Distance)](https://www.r-pkg.org/pkg/Distance)
[![CRAN Version](http://www.r-pkg.org/badges/version/Distance)](https://www.r-pkg.org/pkg/Distance)
 [![Codecov test coverage](https://codecov.io/gh/DistanceDevelopment/Distance/branch/master/graph/badge.svg)](https://codecov.io/gh/DistanceDevelopment/Distance?branch=master)

`Distance` is a simple way of fitting detection functions to distance sampling data for both line and point transects. Adjustment term selection, left and right truncation as well as monotonicity constraints and binning are supported. Abundance and density estimates can also be calculated (via a Horvitz-Thompson-like estimator) if survey area information is provided.

# Using `Distance`

For more information and examples of use [take a look at this paper](https://www.jstatsoft.org/article/view/v089i01) published in Journal of Statistical Software in May 2019.

We also maintain a set of example analyses at [examples.distancesampling.org](http://examples.distancesampling.org).

# Getting `Distance`

The easiest way to ensure you have the latest version of `Distance`, is to install Hadley Wickham's `devtools` package:

      install.packages("devtools")

then install `Distance` from github:

      library(devtools)
      install_github("DistanceDevelopment/Distance")


