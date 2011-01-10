\name{gof}
\alias{gof}

\title{
gof
}
\description{
Performs the nonparametric goodness of fit test on a single
observation. It returns a test statistic for the observation, a vector
of simulated test statistics, and a corresponding p value for a two tailed
test.
}
\usage{
gof(observed, simulated, weights = NULL)
}

\arguments{
  \item{observed}{
A vector valued numeric corresponding to the single observation over which
goodness of fit is to be tested
}
  \item{simulated}{
A \code{gof.preprocess} object returned by \link{gof.preprocess} corresponding to
the null distribution.
}
  \item{weights}{
If null, a default even weighting across the test statistics is used. Otherwise,
a definite square matrix of length equal to the coordinates in the observation
may be supplied. See \link{gof.optimize} for a function to optimize power the 
over the weighting matrix if a particular alternative is in mind.
}
}
\details{
See the working paper by Lospinoso and Snijders (2011) for more information.
}
\value{
list(p=p, v.obs=testStatistic.observed, v.sim=testStatistic.simulated)
\item{p }{p-value corresponding to the probability that, in a hypothetical long run of draws
from the null distribution, we would observe evidence stronger in favor of rejection than
the evidence we have currently observed.}
\item{v.obs }{The test statistic of the observation.}
\item{v.sim }{The test statistics of the simulated draws from the null distribution.}
}
\references{
\url{http://stats.ox.ac.uk/~lospinos}
}
\author{
Josh Lospinoso
}
\seealso{
\link{snopgof}
}
\examples{
# Standard MVN null distribution
null.obs = rnorm(3)

# Simulate 25 draws from the null distribution
null.sim = matrix(rnorm(75),25,3)

# Preprocess the simulations
(null.preproc = gof.preprocess(null.sim))

# Conduct the test
gof(null.obs, null.preproc)
}

\keyword{ htest
nonparametric }

