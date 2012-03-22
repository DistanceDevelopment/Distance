\name{ds}
\alias{ds}
\title{Fit detection functions to line or point transect data}
\usage{
  ds(data, truncation = NULL, transect = "line",
    formula = ~1, key = "hn", adjustment = "cos",
    order = NULL, scale = "scale", cutpoints = NULL,
    monotonicity = FALSE, region.table = NULL,
    sample.table = NULL, obs.table = NULL,
    convert.units = 1)
}
\arguments{
  \item{data}{a \code{data.frame} containing at least a
  column called \code{distance}.}

  \item{truncation}{either truncation distance (numeric,
  e.g. 5) or percentage (as a string, e.g. "15\%"). Can be
  supplied as a \code{list} with elements \code{left} and
  \code{right} if left truncation is required (e.g.
  \code{list(left=1,right=20)} or
  \code{list(left="1\%",right="15\%")} or even When
  specified as a percentage, the largest \code{right} and
  smallest \code{left} percent distances are discarded.}

  \item{transect}{indicates transect type "line" (default)
  or "point".}

  \item{formula}{formula for the scale parameter. For a CDS
  analysis leave this as its default \code{~1}.}

  \item{key}{key function to use; "hn" gives half-normal
  (default) or "hr" gives hazard-rate.}

  \item{adjustment}{adjustment terms to use; "cos" gives
  cosine (default), "herm" gives Hermite polynomial and
  "poly" gives simple polynomial.  "cos" is recommended.}

  \item{order}{orders of the adjustment terms to fit (as a
  vector/scalar), the default value (\code{NULL}) will
  select via AIC. For cosine adjustments, valid orders are
  integers greater than 2. For Hermite polynomials, even
  integers equal or greater than 4 are allowed. For simple
  polynomials even integers equal or greater than 2 are
  allowed.}

  \item{scale}{the scale by which the distances in the
  adjustment terms are divided. Defaults to "scale", the
  scale parameter of the detection function. Other option
  is "width" which scales by the truncation distance.}

  \item{cutpoints}{if the data are binned, this vector
  gives the cutpoints of the bins. Ensure that the first
  element is 0 (or the left truncation distance) and the
  last is the distance to the end of the furthest bin.
  (Default \code{NULL}, no binning.)  Note that if
  \code{data} has columns \code{distbegin} and
  \code{distend} then these will be used as bins if
  \code{cutpoints} is not specified. If both are specified,
  \code{cutpoints} has precedence.}

  \item{monotonicity}{should the detection function be
  constrained for monotonicity weakly ("weak"), strictly
  ("strict") or not at all ("none" or \code{FALSE}). See
  Montonicity, below. (Default \code{FALSE}).}

  \item{region.table}{\code{data.frame} with two columns:
  \tabular{ll}{ \code{Region.Label} \tab label for the
  region\cr \code{Area} \tab area of the region\cr}
  \code{region.table} has one row for each stratum. If
  there is no stratification then \code{region.table} has
  one entry with \code{Area} corresponding to the total
  survey area.}

  \item{sample.table}{\code{data.frame} mapping the regions
  to the samples ( i.e. transects). There are three
  columns: \tabular{ll}{\code{Sample.Label} \tab label for
  the sample\cr \code{Region.Label} \tab label for the
  region that the sample belongs to.\cr \code{Effort} \tab
  the effort expended in that sample (e.g. transect
  length).\cr}}

  \item{obs.table}{\code{data.frame} mapping the individual
  observations (objects) to regions and samples. There
  should be three columns: \tabular{ll}{\code{object} \tab
  \cr \code{Region.Label} \tab label for the region that
  the sample belongs to.\cr \code{Sample.Label} \tab label
  for the sample\cr}}

  \item{convert.units}{conversion between units for
  abundance estimation, see "Units", below. (Defaults to 1,
  implying all of the units are "correct" already.)}
}
\value{
  a list with elements: \tabular{ll}{\code{dsmodel} \tab a
  detection function model object.\cr \code{dht} \tab
  abundance/density information (if survey region data was
  supplied, else \code{NULL}).}
}
\description{
  Some blurb here about distance sampling analysis.
}
\section{Details}{
  If abundance estimates are required the
  \code{data.frame}s \code{region.table},
  \code{sample.table} and \code{obs.table} must be
  supplied.
}

\section{Units}{
  In extrapolating to the entire survey region it is
  important that the unit measurements be consistent or
  converted for consistency.  A conversion factor can be
  specified with the \code{convert.units} variable.  The
  values of \code{Area} in \code{region.table}, must be
  made consistent with the units for \code{Effort} in
  \code{sample.table} and the units of \code{distance} in
  the \code{data.frame} that was analyzed.  It is easiest
  if the units of \code{Area} are the square of the units
  of \code{Effort} and then it is only necessary to convert
  the units of \code{distance} to the units of
  \code{Effort}. For example, if \code{Effort} was entered
  in kilometers and \code{Area} in square kilometers and
  \code{distance} in meters then using
  \code{convert.units=0.001} would convert meters to
  kilometers, density would be expressed in square
  kilometers which would then be consistent with units for
  \code{Area}.  However, they can all be in different units
  as long as the appropriate composite value for
  \code{convert.units} is chosen.  Abundance for a survey
  region can be expressed as: \code{A*N/a} where \code{A}
  is \code{Area} for the survey region, \code{N} is the
  abundance in the covered (sampled) region, and \code{a}
  is the area of the sampled region and is in units of
  \code{Effort * distance}.  The sampled region \code{a} is
  multiplied by \code{convert.units}, so it should be
  chosen such that the result is in the same units as
  \code{Area}.  For example, if \code{Effort} was entered
  in kilometers, \code{Area|} in hectares (100m x 100m) and
  \code{distance} in meters, then using
  \code{convert.units=10} will convert \code{a} to units of
  hectares (100 to convert meters to 100 meters for
  distance and .1 to convert km to 100m units).
}

\section{Monotonicity}{

}
\examples{
# An example from mrds, the golf tee data.
 library(distance)
 data(book.tee.data)
 tee.data<-book.tee.data$book.tee.dataframe
 ds.model<-ds(tee.data,4)
 summary(ds.model)

 # same model, but calculating abundance
 # need to supply the region, sample and observation tables
 region<-book.tee.data$book.tee.region
 samples<-book.tee.data$book.tee.samples
 obs<-book.tee.data$book.tee.obs

 ds.dht.model<-ds(tee.data,4,region.table=region,
              sample.table=samples,obs.table=obs)
 summary(ds.dht.model)

 # specify order 2 cosine adjustments
 ds.model.cos2<-ds(tee.data,4,adjustment="cos",order=2)
 summary(ds.model.cos2)

 # specify order 2 and 4 cosine adjustments
 ds.model.cos24<-ds(tee.data,4,adjustment="cos",order=c(2,4))
 summary(ds.model.cos24)

 # truncate the largest 10\% of the data and fit only a hazard-rate
 # detection function
 ds.model.hr.trunc<-ds(tee.data,truncation="10\%",key="hr",adjustment=NULL)
 summary(ds.model.hr.trunc)
}
\author{
  David L. Miller
}

