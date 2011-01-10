\name{snopgof-package}
\alias{snopgof-package}
\alias{snopgof}
\docType{package}
\title{
snopgof -- Single observation, NOn-Parametric Goodness Of Fit
}
\description{
Implementation of a Single observation, NOn-Parametric Goodness Of Fit 
(SNOPGOF) designed for testing GOF with a single observation supposed 
to come from a distribution with no available closed form 
(but can be easily simulated).
}
\details{
\tabular{ll}{
Package: \tab snopgof\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-01-05\cr
License: \tab GNU GPL\cr
LazyLoad: \tab yes\cr
}

}
\author{
Josh Lospinoso <lospinos@stats.ox.ac.uk>
}
\references{
This package represents a work in progress reflecting a working paper by
Lospinoso and Snijders 2011. Contact the author for details.
}

\keyword{ package 
htest
nonparametric}
\seealso{
\link{gof.compare}
\link{gof.optimize}
\link{gof.preprocess}
\link{gof}
}
\examples{
## In this example, we generate a simulated `null' set. We then preprocess this
## null set. Next, a simulated `alternate' set is used to calibrate a weighting
## matrix by maximizing power to detect this alternate set. Finally, snopgof
## is tested to produce receiver operating characteristic (ROC) curves by
## simulating more draws from the null and alternate distributions and performing
## the test. This file is currently set up to draw from MVN distributions.

## Note that these values are set to low numbers to make passing checks quicker.
## We highly recommend increasing the number of simulations.

## How many coordinates in the MVN distributions?
variates = 3

## Number of simulations to use for null and alternate sets:
simulations.sets = 50

## and for constructing the ROC curves:
simulations.test = 25

## Parameters of the null MVN distribution:
null.mean = runif(variates,-10,10)
temp = matrix(runif(variates^2,-1,1),variates)
null.cov = temp \%*\% t(temp)
null.chol = chol(null.cov)

## Parameters of the null MVN distribution:
alt.mean = runif(variates,0,10)
temp = matrix(runif(variates^2,-1,1),variates)
alt.cov = temp \%*\% t(temp)
alt.chol = chol(alt.cov)

## Generate the null and alternate simulation sets 
null.sims <- matrix(rnorm(variates*simulations.sets), simulations.sets) \%*\% null.chol + matrix(null.mean,simulations.sets,variates,byrow=TRUE)
alt.sims <- matrix(rnorm(variates*simulations.sets), simulations.sets) \%*\% alt.chol + matrix(alt.mean,simulations.sets,variates,byrow=TRUE)

# Pre-process the comparisons for the simulations.
# See ?gof.preprocess for more information.
(null.preproc <- gof.preprocess(null.sims))

# Optimize the power of the test for a simulated alternative
# See ?gof.optimize for more information.
(weights <- gof.optimize(null.preproc, alt.sims))

## Generate the null and alternate observations for the ROC curve
null.obs <- matrix(rnorm(variates*simulations.test), simulations.test) \%*\% null.chol + matrix(null.mean,simulations.test,variates,byrow=TRUE)
alt.obs <- matrix(rnorm(variates*simulations.test), simulations.test) \%*\% alt.chol + matrix(alt.mean,simulations.test,variates,byrow=TRUE)

# Do a tail-type test
# See ?gof for more information. Note that it is wrapped by sapply in order to construct the ROC curve
(res.null <- sapply(1:simulations.test, function(i) gof(null.obs[i,], null.preproc, weights)$p))

## Do a tail-type test by drawing from the null distribution
(res.alt <- sapply(1:simulations.test, function(i) gof(alt.obs[i,], null.preproc, weights)$p))

## Plot the ROC curves
yaxis = seq(0,1,along.with=res.null)
plot(yaxis~sort(res.null), type="l", main="Reciever Operating Characteristic (ROC) Curves", 
	xlab="False Positive", ylab="True Positive (power)")
lines(yaxis~sort(res.alt), type="l", col="red")
legend(.25~.75, c("Null","Alternate"), pch=c(19,19), col=c(1,2))
}
