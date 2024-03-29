% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{generateDiscreteData}
\alias{generateDiscreteData}
\title{A data generation function using a discrete distribution for Rankin
score rather than a normal distribution}
\usage{
generateDiscreteData(prevalence, N, support = 0L:6L, ctlDist, trtDist)
}
\arguments{
\item{prevalence}{a vector of group prevalences (length denoted by J below)}

\item{N}{the sample size to generate}

\item{support}{the support values of the discrete distribution (length K), default 0:6}

\item{ctlDist}{a probability vector of length K denoting the Rankin score distribution for control.}

\item{trtDist}{an K x J probability matrix with each column is the Rankin distribution for the associated group}
}
\value{
a three-column data frame of \code{subGroup}, \code{trt} (0 or 1), and \code{score}
}
\description{
A data generation function using a discrete distribution for Rankin
score rather than a normal distribution
}
\examples{
# Simulate data from a discrete distribution for the Rankin scores,
# which are typically ordinal integers from 0 to 6 in the following
# simulations. So we define a few scenarios.
library(ASSISTant)
null.uniform <- rep(1, 7L) ## uniform on 7 support points
hourglass <- c(1, 2, 2, 1, 2, 2, 1)
inverted.hourglass <- c(2, 1, 1, 2, 1, 1, 2)
bottom.heavy <- c(2, 2, 2, 1, 1, 1, 1)
bottom.heavier <- c(3, 3, 2, 2, 1, 1, 1)
top.heavy <- c(1, 1, 1, 1, 2, 2, 2)
top.heavier <- c(1, 1, 1, 2, 2, 3, 3)
ctlDist <- null.uniform
trtDist <- cbind(null.uniform, null.uniform, hourglass, hourglass) ## 4 groups
generateDiscreteData(prevalence = rep(1, 4), N = 10, ctlDist = ctlDist,
                     trtDist = trtDist) ## default support is 0:6
trtDist <- cbind(bottom.heavy, bottom.heavy, top.heavy, top.heavy)
generateDiscreteData(prevalence = rep(1, 4), N = 10, ctlDist = ctlDist,
                     trtDist = trtDist)
support <- c(-2, -1, 0, 1, 2) ## Support of distribution
top.loaded <- c(1, 1, 1, 3, 3) ## Top is heavier
ctl.dist <- c(1, 1, 1, 1, 1) ## null on 5 support points
trt.dist <- cbind(ctl.dist, ctl.dist, top.loaded) ## 3 groups
generateDiscreteData(prevalence = rep(1, 3), N = 10, support = support,
                     ctlDist = ctl.dist, trtDist = trt.dist)
## ctl.dist can also be a matrix with different nulls for each subgroup
uniform <- rep(1, 5)
bot.loaded <- c(3, 3, 1, 1, 1)
ctl.dist <- matrix(c(uniform, bot.loaded, top.loaded), nrow = 5)
generateDiscreteData(prevalence = rep(1, 3), N = 10, support = support,
                     ctlDist = ctl.dist, trtDist = trt.dist)
}
