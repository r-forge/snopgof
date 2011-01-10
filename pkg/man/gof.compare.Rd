\name{gof.compare}
\alias{gof.compare}

\title{
gof.compare
}
\description{
Used internally to \link{snopgof} to perform the comparisons which form the
tail statistics.
}
\usage{
gof.compare(observed, simulated, centralityMeasure="mean")
}

\arguments{
  \item{observed}{
Vector valued numeric corresponding to the observation.
Must have the same length as simulated has columns.
}
  \item{simulated}{
Simulated draws from the null distribution. Each row corresponds to
a draw.
}
  \item{centralityMeasure}{
Which centrality measure to use. Accepts "mean", "median", and "mode."
}
}
\details{
This function should not need to be called by the typical user.
}
\value{
A matrix of the comparisons for the tail statistics. Each row corresponds to a tail statistic
which can be interpreted, for row i, as "How many draws from the null observation had values along exactly i 
coordinates which were closer to the chosen centrality measure (i.e. mean, median, or mode) than this observation?" Thus,
each column corresponds to an observation.
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
# See the source of \link{gof.preprocess} for advanced uses of this function.
}

\keyword{ htest
nonparametric }
