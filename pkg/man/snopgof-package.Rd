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
\link{gof.mhd.asymptotic}
\link{gof}
}
\examples{
\dontrun{
## In this example, we generate a simulated `null' set. We then preprocess this
## null set. Next, a simulated `alternate' set is used to calibrate a weighting
## matrix by maximizing power to detect this alternate set. Finally, snopgof
## is tested to produce receiver operating characteristic (ROC) curves by
## simulating more draws from the null and alternate distributions and performing
## the test.

## The file is set up to draw independent Cauchy random variables with a large
## location shift in the alternate distribution. These are chosen to illustrate
## a case where the proposed test of Lospinoso and Snijders (2011) outperforms
## the Mahalanobis distance test statistic.

## How many coordinates to include:
variates = 3
## Number of simulations to use for null and alternate sets:
simulations.sets = 1500

## Generate the null and alternate simulation sets 
null.sims <- matrix(rcauchy(simulations.sets*variates), ncol=variates)
alt.sims <- matrix(rcauchy(simulations.sets*variates, location=10, scale=1), ncol=variates)

## Generate observations from the same distributions
null.obs <- matrix(rcauchy(simulations.sets*variates), ncol=variates)
alt.obs <- matrix(rcauchy(simulations.sets*variates, location=10, scale=1), ncol=variates)

# Pre-process the comparisons for the simulations.
# See ?gof.preprocess for more information.
(null.preproc <- gof.preprocess(null.sims))

# Optimize the power of the test for a simulated alternative
# See ?gof.optimize for more information.
(weights <- gof.optimize(null.preproc, alt.sims, sannIterations=50, verbose=TRUE))

# Do a tail-type test while drawing from the null distribution
# See ?gof for more information. Note that it is wrapped by sapply in order to construct the ROC curve
res.null <- gof(null.obs, null.preproc, weights)

## Do a tail-type test while drawing from the alternate distribution
res.alt <- gof(alt.obs, null.preproc, weights)

## Plot the ROC curves
yaxis = seq(0,1,along.with=res.null$p)
plot(yaxis~sort(res.null$p), type="l", main="Receiver Operating Characteristic (ROC) Curves", 
	xlab="False Positive", ylab="True Positive (power)")

lines(yaxis~sort(res.null$p.MHD), type="l", col=2)
lines(yaxis~sort(res.alt$p), type="l", col=3)
lines(yaxis~sort(res.alt$p.MHD), type="l", col=4)
legend(.25~.75, c("Null","Null (Mah. Dist.", "Alt.","Alt. (Mah. Dist.)"), pch=c(19,19,19,19), col=c(1,2,3,4))
}

#################################################################################
## Instead of constructing ROC curves, do a simple single observation test:
## How many coordinates in the MVN distributions?
variates = 10
simulations.sets = 500

## Parameters of the null MVN distribution:
null.mean = runif(variates,-10,10)
temp = matrix(runif(variates^2,-1,1),variates)
null.cov = temp \%*\% t(temp)
null.chol = chol(null.cov)

## Generate the null and alternate simulation sets
null.sims <- matrix(rnorm(variates*simulations.sets), simulations.sets) \%*\% null.chol + matrix(null.mean,simulations.sets,variates,byrow=TRUE)
null.obs <- as.numeric(rnorm(variates) \%*\% null.chol + null.mean)
pre <- gof.preprocess(null.sims)
(results <- gof(null.obs, pre))
plot(results)
}
