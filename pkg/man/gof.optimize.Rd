\name{gof.optimize}
\alias{gof.optimize}

\title{
gof.optimize
}
\description{
An optional function which can be used to optimize the weighting matrix against some simulated
alternative for a specified type one error level.
}

\usage{
gof.optimize(null, alternate, typeOneError = 0.05, weights = NULL, 
	verbose=FALSE, sannIterations=1000, permIterations=1023)
}

\arguments{
  \item{null}{
  A \cite{gof.preprocess} object containing the preprocessed simulated draws from the null distribution.
}
  \item{alternate}{
 A matrix of draws from the desired alternate distribution. Power will be maximized to detect draws from
 this distribution at the specified type one error level. The coordinates correspond to columns so that
 each draw corresponds to a row.
}
  \item{typeOneError}{
 The desired level of type one error, i.e., at which level will you decide to reject the null hypothesis
 in favor of the alternative.
}
  \item{weights}{
 An initial value for the optimization problem.
}
  \item{verbose}{
 Whether to report more or less incremental information while the annealing is optimizing.
}
 \item{sannIterations}{
 How many iterations of simulated annealing to attempt. If zero, only the permutations will be considered
}
\item{permIterations}{
How many permutations to consider, at a maximum. Meaningful thresholds are necessarily in the set 1 3 7 15 
31 63 127 255 511 1023 2047 4095 8191 16383 32767 65535 131071 262143 524287 1048575. The default is set
to exhaust ten coordinates. When the threshold is hit, the most inclusive permutations are included first.
}
}
\details{
The \link{optim} function is used with the simulated annealing option ("SANN") to optimize power for the
given type one error. The results may or may not be satistfactory, and may require some tweaking.
This function will try all combinations of diagonal elements as weights in an attempt to find a good
solution. 
}

\value{
A weighting matrix which may be used as input to \link{gof}
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
# See ?snopgof 
}

\keyword{ htest
nonparametric }
