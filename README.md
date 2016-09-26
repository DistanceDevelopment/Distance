`Distance`
==========

[![Build Status](https://travis-ci.org/DistanceDevelopment/Distance.svg?branch=master)](https://travis-ci.org/DistanceDevelopment/Distance)
[![CRAN (RStudio Mirror) Downloads](http://cranlogs.r-pkg.org/badges/Distance)](http://www.r-pkg.org/pkg/Distance)
[![CRAN Version](http://www.r-pkg.org/badges/version/Distance)](http://www.r-pkg.org/pkg/Distance)

`Distance` is a simple way of fitting detection functions to distance sampling data for both line and point transects. Adjustment term selection, left and right truncation as well as monotonicity constraints and binning are supported. Abundance and density estimates can also be calculated (via a Horvitz-Thompson-like estimator) if survey area information is provided.

# Using `Distance`

For more information and examples of use [take a look at this preprint](http://biorxiv.org/content/early/2016/07/14/063891).

# Getting `Distance`

The easiest way to ensure you have the latest version of `Distance`, is to install Hadley Wickham's devtools package:

      install.packages("devtools")

then install `Distance` from github:

      library(devtools)
      install_github("DistanceDevelopment/Distance")


