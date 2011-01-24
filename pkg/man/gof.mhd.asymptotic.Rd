\name{gof.mhd.asymptotic}
\alias{gof.mhd.asymptotic}

\title{
gof.mhd.asymptotic
}
\description{
Performs a Mahalanobis distance type goodness of fit test on a single
observation. It returns a test statistic for the observation and a 
corresponding p value.
}
\usage{
gof.mhd.asymptotic(observed, null.mean, null.cov, twoTailed = FALSE)
}

\arguments{
  \item{observed}{
A vector valued numeric corresponding to the single observation over which
goodness of fit is to be tested
}
  \item{null.mean}{
A vector valued numeric corresponding to the mean of the MVN null distribution.
}
  \item{null.cov}{
A matrix value corresponding to the covariance of the MVN null distribution.
}
  \item{twoTailed}{
A boolean indicating whether to conduct a two-tailed test on the Chi-squared distribution.
}
}
\details{
See the working paper by Lospinoso and Snijders (2011) for more information.
}
\value{
list(p=p, v.obs=testStatistic.observed)
\item{p }{p-value corresponding to the probability that, in a hypothetical long run of draws
from the null distribution, we would observe evidence stronger in favor of rejection than
the evidence we have currently observed.}
\item{v.obs }{The test statistic of the observation.}
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
# Standard MVN null distribution, Receiver Operating Characteristic curve
null.mean = rep(0, 5)
null.cov = matrix(0,5,5)
diag(null.cov) <- 1
x<-sapply(1:5000, function (i) {gof.mhd.asymptotic(rnorm(5), null.mean, null.cov)$p})
plot(sort(x),main="ROC Curve for null distributed data")
}

\keyword{ htest
nonparametric }

