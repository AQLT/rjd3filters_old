% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/x11.R
\name{x11}
\alias{x11}
\title{X-11 Decomposition With Custom Trend Filters}
\usage{
x11(
  x,
  period = frequency(x),
  trend.coefs,
  mul = TRUE,
  seas.s0 = c("S3X3", "S3X1", "S3X5", "S3X9", "S3X15"),
  seas.s1 = c("S3X5", "S3X3", "S3X1", "S3X9", "S3X15"),
  extreme.lsig = 1.5,
  extreme.usig = 2.5,
  userdefined = NULL
)
}
\arguments{
\item{x}{input time-series.}

\item{period}{period.}

\item{trend.coefs}{coefficients of the filters used for the trend-cycle extraction from
the real-time asymmetric filter to the symmetric filter. Can be a, object of class \code{"list"},
\code{"matrix"}, \code{"lp_filter"} or \code{"rkhs_filter"}.}

\item{mul}{boolean indicating if the decomposition mode is multiplicative.}

\item{seas.s0, seas.s1}{seasonal filters.}

\item{extreme.lsig, extreme.usig}{boundaries used for outlier correction in irregular.}

\item{userdefined}{a vector containing the additional output variables.}
}
\description{
Perform the X-11 decomposition using custom trend filter
}
\examples{
x <- retailsa$AllOtherGenMerchandiseStores
decomposition_lp <- x11(x, trend.coefs = lp_filter())
decomposition_rkhs <- x11(x, trend.coefs = rkhs_filter())
plot(x)
lines(decomposition_lp$decomposition[,"t"], col = "red")
lines(decomposition_rkhs$decomposition[,"t"], col = "green")
}
