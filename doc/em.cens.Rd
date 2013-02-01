\name{em.cens}
\alias{em.cens}
\title{An EM implementation for fitting censored mixtures with maximum likelihood}
\description{
  attempts to fit a k-component censored mixture using EM
}
\usage{ em.cens(y, i0.mat, w0, m0, sd0, niter=30)
}
\arguments{
  \item{y}{a n-length numeric vector containing the dataset}
    \item{i0.mat}{a n by k matrix of 1/0 flags indicating an initial guess for which component each observation belongs to}
    \item{w0}{a k-length numeric vector of initial component weights}
    \item{m0}{a k-length numeric vector of initial component positions}
    \item{sd0}{a k-length numeric vector of initial component standard deviations}
    \item{niter}{number of iterations of EM to run}
}
\value{
    list(id=, wt=, mean=, sd=)
    \item{id}{a n by k matrix of probability that observation m belongs to component j}
    \item{wt}{a k-length numeric vector of current component weight estimates}
    \item{mean}{a k-length numeric vector of current component position estimates}
    \item{sd}{a k-length numeric vector of current component standard deviation estimates}
    
}
