\name{Sbi.ureda}
\alias{Sbi.ureda}
\title{a weighted version of a robust measure of scale, Sbi, as described in "Understanding robust and exploratory data analysis"}
\description{
  A weighted version of Sbi, particularly useful for mixture models
}
\usage{ Sbi.ureda(x, bi, c0, wt, s0=0, verbose=FALSE)
}
\arguments{
  \item{x}{a n-length vector of data}
  \item{bi}{the location point around which the scale is measured (bi because Sbi is associated with the biweight)}
  \item{c0}{a tuning parameter}
  \item{wt}{a n-length vector of weights}
  \item{s0}{a starting value for the scale, defaults to a weighted median absolute deviation from bi}
  \item{verbose}{if TRUE, will print many of the intermediate steps of the calculation}

}

\references{
Understanding robust and exploratory data analysis
edited by: David C. Hoaglin, Frederick Mosteller, John W. Tukey
(1983) 
}

