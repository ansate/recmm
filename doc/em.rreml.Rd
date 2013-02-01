\name{em.rreml}
\alias{em.rreml}
\title{A simple EM implementation for finding the parameters of Gaussian mixtures, in the case of censored data}
\description{
  attempts to fit a k-component mixture using EM, taking into account censoring, and having the option to use robust measures
}
\usage{ em.rreml(y, i0.mat, w0, m0, sd0, niter=30, verbose=FALSE, robust=FALSE)
}
\arguments{
  \item{y}{a n-length numeric vector containing the dataset}
    \item{i0.mat}{a n by k matrix of 1/0 flags indicating an initial guess for which component each observation belongs to}
    \item{w0}{a k-length numeric vector of initial component weights}
    \item{m0}{a k-length numeric vector of initial component positions}
    \item{sd0}{a k-length numeric vector of initial component standard deviations}
    \item{niter}{number of iterations of EM to run}
    \item{verbose}{print details of the iterations}
    \item{robust}{use the biweight and Sbi, robust measures of location and spread}
}
\value{
    list(id=, wt=, mean=, sd=)
    \item{id}{a n by k matrix of probability that observation m belongs to component j}
    \item{wt}{a k-length numeric vector of current component weight estimates}
    \item{mean}{a k-length numeric vector of current component position estimates}
    \item{sd}{a k-length numeric vector of current component standard deviation estimates}
    
}
