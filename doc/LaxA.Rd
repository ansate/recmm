\name{LaxA}
\alias{LaxA}
\title{a weighted version of a robust measure of scale, A, as described in "Robust Estimators of Scale: Finite-Sample Performance in Long-Tailed Symmetric Distributions"}
\description{
  A weighted version of Lax's A, particularly useful for mixture models
}
\usage{ LaxA(x, bi, c0, wt, s0=0)
}
\arguments{
  \item{x}{a n-length vector of data}
  \item{bi}{the location point around which the scale is measured (bi because Sbi is associated with the biweight)}
  \item{c0}{a tuning parameter}
  \item{wt}{a n-length vector of weights}
  \item{s0}{a starting value for the scale, defaults to a weighted median absolute deviation from bi}

}

\references{
Robust Estimators of Scale: Finite-Sample Performance in Long-Tailed Symmetric Distributions
by: David A. Lax
Journal of the American Statistical Association, Vol. 80, No. 391. (1985), pp. 736-741.
}

