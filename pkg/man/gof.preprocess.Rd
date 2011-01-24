\name{gof.preprocess}
\alias{gof.preprocess}
\title{
gof.preprocess
}
\description{
Preprocesses the comparisons within a simulated set of draws from the
null distribution. These comparisons are used to construct a simulated
distribution of tail statistics. Since this needs to be done only once,
it is abstracted from \link{gof.optimize} and \link{gof} for computational
reasons.
}
\usage{
gof.preprocess(simulated, creasingFunction=median)
}

\arguments{
  \item{simulated}{
A matrix of simulated draws from the null distribution.
Each column should represent a coordinate so that each row
corresponds to an independent draw from the null distribution.
}
\item{creasingFunction}{Which function is used to center the simulations.}
}

\value{
\item{Comparisons }{A matrix of the comparisons for the tail statistics. Each row corresponds to a tail statistic
which can be interpreted, for row i, as "How many draws from the null observation had values along exactly i 
coordinates which were closer to the chosen centrality measure (i.e. mean, median, or mode) than this observation?" Thus,
each column corresponds to an observation.}
\item{Iterations }{The number of iterations comprising the simulated draws from the null distribution.}
\item{Variates }{The number of coordinates in the support of the null distribution.}
\item{Original }{The original simulations given from the null distribution, i.e. \code{simulated}.}
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
# See ?snopgof for a full description.

# Standard MVN null distribution
m=1000; v=5;
null.obs = rnorm(v)

# Simulate 25 draws from the null distribution
null.sim = matrix(rnorm(m*v),m,v)

# Preprocess the simulations
(null.preproc = gof.preprocess(null.sim))

# Conduct the test
(x<-gof(null.obs, null.preproc))

# Plot the results
plot(x)
}

\keyword{ htest
nonparametric }
