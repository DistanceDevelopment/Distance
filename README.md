Distance: analysis of distance sampling data

[![R-CMD-check](https://github.com/DistanceDevelopment/Distance/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/DistanceDevelopment/Distance/actions/workflows/check-standard.yaml)
[![CRAN (RStudio Mirror) Downloads](http://cranlogs.r-pkg.org/badges/Distance)](https://www.r-pkg.org/pkg/Distance)
[![CRAN Version](http://www.r-pkg.org/badges/version/Distance)](https://www.r-pkg.org/pkg/Distance)
 [![Codecov test coverage](https://codecov.io/gh/DistanceDevelopment/Distance/branch/master/graph/badge.svg)](https://app.codecov.io/gh/DistanceDevelopment/Distance?branch=master)

`Distance` is a simple way of fitting detection functions to distance sampling data for both line and point transects. Adjustment term selection, left and right truncation as well as monotonicity constraints and binning are supported. Abundance and density estimates can also be calculated (via a Horvitz-Thompson-like estimator) if survey area information is provided.

# Using `Distance`

### Distance **R** package preferred citation
- Miller, D. L., Rexstad, E., Thomas, L., Marshall, L., & Laake, J. L. (2019). Distance Sampling in R. Journal of Statistical Software, 89(1), 1–28. DOI: [10.18637/jss.v089.i01](https://doi.org/10.18637/jss.v089.i01)

Consult the [Articles](https://distancedevelopment.github.io/Distance/articles/) section of this site for case studies of distance sampling analyses.

# Getting `Distance`

The easiest way to ensure you have the latest version of `Distance`, is to install `devtools`:

```{r}
install.packages("devtools")
```

then install `Distance` from Github:

```{r}
library(devtools)
install_github("DistanceDevelopment/Distance")
```