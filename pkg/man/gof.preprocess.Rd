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
gof.preprocess(simulated)
}

\arguments{
  \item{simulated}{
A matrix of simulated draws from the null distribution.
Each column should represent a coordinate so that each row
corresponds to an independent draw from the null distribution.
}
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
# Standard MVN null distribution
null.obs = rnorm(3)

# Simulate 25 draws from the null distribution
null.sim = matrix(rnorm(75),25,3)

# Preprocess the simulations
(null.preproc = gof.preprocess(null.sim))

}

\keyword{ htest
nonparametric }
