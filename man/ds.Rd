\name{ds}
\alias{ds}
\title{Fit detection functions to line or point transect data}
\usage{
  ds(data, truncation, transect = "line", formula = ~1,
    key = "hn", adjustment = "cos", order = NULL,
    scale = "scale", binned = FALSE)
}
\arguments{
  \item{data}{a \code{data.frame} containing at least a
  column called \code{distance}.}

  \item{truncation}{truncation distance.}

  \item{transect}{indicates transect type "line" (default)
  or "point".}

  \item{formula}{formula for the scale parameter. For a CDS
  analysis leave this as its default \code{~1}.}

  \item{key}{key function to use; "hn" gives half-normal
  (default) or "hr" gives hazard-rate.}

  \item{adjustment}{adjustment terms to use; "cos" gives
  cosine (default), "herm" gives Hermite polynomial and
  "poly" gives simple polynomial.  "cos" is recommended.}

  \item{order}{orders of the adjustment terms (as a
  vector), the default value (\code{NULL}) will select via
  AIC. TKTKTK talk about which orders are correct}

  \item{scale}{the scale by which the distances in the
  adjustment terms are divided. Defaults to "scale", the
  scale parameter of the detection function. Other option
  is "width" which scales by the truncation distance.}

  \item{binned}{if the data are binned, set this to
  \code{TRUE} (default \code{FALSE}).}
}
\description{
  BLURB
}
\examples{
# An example from mrds, the golf tee data.
 library(distance)
 data(book.tee.data)
 tee.data<-book.tee.data$book.tee.dataframe
 ds.model<-ds(tee.data,4)
 summary(ds.model)
}
\author{
  David L. Miller
}

