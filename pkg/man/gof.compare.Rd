\name{gof.compare}
\alias{gof.compare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gof.compare
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gof.compare(observed, simulated, median = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{observed}{
%%     ~~Describe \code{observed} here~~
}
  \item{simulated}{
%%     ~~Describe \code{simulated} here~~
}
  \item{median}{
%%     ~~Describe \code{median} here~~
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
function(observed, simulated, median=TRUE) {
if (class(observed) != "numeric" | class(simulated) != "matrix") {
stop("Parameters are invalid.")
  }
if (length(observed) != dim(simulated)[2]) {
stop("Dimensionality of function parameters do not match.")
  }
variates = length(observed)
simulations = dim(simulated)[1]
if (median) {
centroid = apply(simulated,2,median)
  } else {
centroid = apply(simulated,2,mean)
  }
sapply(0:(variates-1), function (threshold) {sum(sapply(1:simulations, 
function (i) sum(abs(observed - centroid) > 
abs(simulated[i,] - centroid)) == threshold)) 
  })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
