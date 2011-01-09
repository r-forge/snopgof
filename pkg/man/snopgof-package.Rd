\name{snopgof-package}
\alias{snopgof-package}
\alias{snopgof}
\docType{package}
\title{
What the package does (short line)
~~ package title ~~
}
\description{
More about what it does (maybe more than one line)
~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab snopgof\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-01-05\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
library(snopgof)
## Length of an observation
variates = 5
## Number of simulations to use for EDF
simulations.null = 500
simulations.alt = 250
simulations.test = 250
## Alternate normal distribution, mean
alt.shift = 0
## Alternate normal distribution, variance
alt.scale = .75

## Generate an observation of length variates ~ Normal(0,I) where I: Identity matrix
obs <- rnorm(variates)
## Generate 'simulations' observations of length variates ~ Normal(0,I) where I: Identity matrix
sim.null <- matrix(rnorm(simulations.null * variates), simulations.null,  variates)
sim.alt <- matrix(rnorm(simulations.alt * variates, alt.shift, alt.scale), simulations.alt,  variates)

# Pre-process the comparisons for the simulations
(pre.null <- gof.preprocess(sim.null))

# Optimize the power of the test for a simulated alternative
(weights <- gof.optimize(pre.null, sim.alt))

## Do a tail-type test when null true
(res.null <- sapply(1:simulations.test, function(i) gof(rnorm(variates), pre.null, weights)$p))

## Do a tail-type test when null false
(res.alt <- sapply(1:simulations.test, function(i) gof(rnorm(variates, alt.shift, alt.scale), pre.null, weights)$p))

## Plot EDF values
#pdf("~/Desktop/Results3.pdf")
plot(sort(res.null), type="l", main="EDF of p-values for null simulations")
lines(sort(res.alt), type="l", col="red")
#dev.off()
}
