\name{gof}
\alias{gof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gof
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gof(observed, simulated, weights = rep(1/null$Variates, null$Variates))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{observed}{
%%     ~~Describe \code{observed} here~~
}
  \item{simulated}{
%%     ~~Describe \code{simulated} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(observed, simulated, weights=rep(1/null$Variates, null$Variates)) {
if (class(observed) != "numeric" | class(simulated) != "gof.preprocess") {
stop("Parameters are invalid.")
  }
if (length(observed) != simulated$Variates) {
stop("Dimensionality of function parameters do not match.")
  }
tailComparisons.observed = gof.EDFTest.compare(observed, simulated$Original)
testStatistic.observed = as.numeric(tailComparisons.observed \%*\% weights)
testStatistic.simulated = as.numeric(weights \%*\% simulated$Comparisons)
edf.observed = sum(testStatistic.simulated <= testStatistic.observed) / simulated$Iterations
p <- 1 - abs(1 - 2 * edf.observed)
return = list(p=p, v.obs=testStatistic.observed, v.sim=testStatistic.simulated)
class(return) <- "gof.EDFTestResult"
return
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
