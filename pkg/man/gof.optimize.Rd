\name{gof.optimize}
\alias{gof.optimize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gof.optimize
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gof.optimize(null, alternate, typeOneError = 0.05, weights = rep(1/null$Variates, null$Variates))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{null}{
%%     ~~Describe \code{null} here~~
}
  \item{alternate}{
%%     ~~Describe \code{alternate} here~~
}
  \item{typeOneError}{
%%     ~~Describe \code{typeOneError} here~~
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
function(null, alternate, typeOneError=.05, weights=rep(1/null$Variates, null$Variates)) {
if (class(null) != "gof.preprocess" | class(alternate) != "matrix" | typeOneError <= 0 | typeOneError >= 1) {
stop("Parameters are invalid.")
  }
if (null$Variates != ncol(alternate) | null$Variates != length(weights)) {
stop("Variates do not match in null and alternate")
  }
fn <- function(w) {
w <- w / sum(w)
a <- sapply(1:nrow(alternate), function(i) gof.EDFTest(alternate[i,],null,w)$p)
cat("Testing power, summary of p-values for weighting ",w,":\n ")
power = sum(a <= typeOneError) / nrow(alternate)
cat(" w. alpha .01 = ", sum(a < .01) / nrow(alternate), "\n",
" w. alpha .05 = ", sum(a < .05) / nrow(alternate), "\n",
" w. alpha .10 = ", sum(a < .10) / nrow(alternate), "\n",
" w. alpha .25 = ", sum(a < .25) / nrow(alternate), "\n")
power
  }
lower = rep(0.000001, null$Variates)
upper = rep(1, null$Variates)
r <- optim(weights,fn,lower=lower,upper=upper,method="L-BFGS-B")
(weights <- r$par / sum(r$par))
weights[weights < .001] = 0
(weights <- r$par / sum(r$par))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
