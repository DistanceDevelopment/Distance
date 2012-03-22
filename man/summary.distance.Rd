\name{summary.distance}
\alias{summary.distance}
\title{Summary of distance sampling analysis}
\usage{
  \method{summary}{distance} (object, ...)
}
\arguments{
  \item{object}{a distance analysis}

  \item{\dots}{unspecified and unused arguments for S3
  consistency}
}
\value{
  list of extracted and summarized objects
}
\description{
  Provides a brief summary of a distance sampling analysis.
  This includes
}
\note{
  This function just calls \code{\link{summary.ds}} and
  \code{\link{dht}}, collates and prints the results in a
  nice way.
}
\author{
  David L. Miller
}
\keyword{utility}

