% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsdata.R
\docType{data}
\name{minke}
\alias{minke}
\title{Simulated minke whale data}
\format{
\code{data.frame} with 99 observations of 5 variables:
\itemize{
\item \code{Region.Label} stratum label (\code{"North"} or \code{"South"})
\item \code{Area} stratum area
\item \code{Sample.Label} transect identifier
\item \code{Effort} transect length
\item \code{distance} observed distance
\item \code{object} unique object ID
}
}
\source{
Shipped with the Distance for Windows.
}
\description{
Data simulated from models fitted to 1992/1993 Southern Hemisphere minke
whale data collected by the International Whaling Commission. See Branch and
Butterworth (2001) for survey details (survey design is shown in figure
1(e)). Data simulated by David Borchers.
}
\details{
Data are included here as both R data and as an Excel spreadsheet to
illustrate the "flat file" input method. See \code{\link{flatfile}} for how
to load this data and an example analysis.
}
\examples{
data(minke)
head(minke)
}
\references{
Branch, T.A. and D.S. Butterworth (2001) Southern Hemisphere
minke whales: standardised abundance estimates from the 1978/79 to 1997/98
IDCR-SOWER surveys. Journal of Cetacean Research and Management 3(2):
143-174

Hedley, S.L., and S.T. Buckland. Spatial Models for Line Transect Sampling.
Journal of Agricultural, Biological, and Environmental Statistics 9, no. 2
(2004): 181-199. \doi{10.1198/1085711043578}.
}
\keyword{datasets}
