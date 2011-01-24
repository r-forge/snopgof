gof <-
function(observed, simulated, weights=NULL) {
	if (class(observed) == "numeric") {
		observed <- matrix(observed,nrow=1)
	}

	if (class(observed) != "matrix" | class(simulated) != "gof.preprocess") {
		stop("simulated must be a gof.preprocess object.")
	}
	if (ncol(observed) != simulated$Variates) {
		stop("Dimensionality of function parameters do not match.")
	}
	observations = nrow(observed)
	
	if (is.null(weights)) {
		if (length(observed) > 1) {
			weights = matrix(0, simulated$Variates+1, simulated$Variates+1)
			diag(weights) <- 1 / (simulated$Variates+1)
		} else {
			weights = 1
		}
	}
	
	foldedObservations <- sapply(1:observations, function(i) abs(observed[i,]-simulated$Expectation))
	p <- .Call("doTest", foldedObservations, simulated$FoldedSimulations, simulated$CenteredComparisons, simulated$ComparisonCrease, weights)
	
	#	Comparisons in R instead of C
	#	tailComparisons.observed = sapply(1:observations, function (i) gof.compare(abs(observed[i,]-simulated$Expectation), simulated$FoldedSimulations) )
	#	tailComparisons.observed.centered = abs(tailComparisons.observed - simulated$ComparisonCrease)
	#	testStatistic.observed =  as.numeric(sapply(1:observations, function(i) t(tailComparisons.observed.centered[,i]) %*% weights %*% tailComparisons.observed.centered[,i]))
	#	testStatistic.simulated = as.numeric(sapply(1:simulated$Iterations, function(i) t(simulated$CenteredComparisons[,i]) %*% weights %*% simulated$CenteredComparisons[,i]))
	#	edf.observed = sapply(1:observations, function (i) sum(testStatistic.simulated <= testStatistic.observed[i]) / simulated$Iterations )
	#	p <- 1 - abs(1 - 2 * edf.observed)
	mhd <- sapply(1:observations, function(i) t(observed[i,] - simulated$Expectation) %*% simulated$InverseCovariance %*% (observed[i,] - simulated$Expectation) )
	p.mhd <- sapply(1:observations, function (i) 1 - abs(1 - 2 * sum(mhd[i] >= simulated$MHDistances)/simulated$Iterations) )
	return = list(p=p, #v.obs=testStatistic.observed, v.sim=testStatistic.simulated,
			Observation=observed, Preprocess=simulated, MHDistance=mhd, p.MHD = p.mhd)
	class(return) <- "gof.result"
	return
}

print.gof.result <- function (x, ...) {
	with(x, cat("snopgof results\nVariates:", Preprocess$Variates ,"\nSimulations:", Preprocess$Iterations ,"\nP-Value:", p,"\nP-Value(MHD):", p.MHD))
}

plot.gof.result <- function (x, ...) {
	boxplot(as.numeric(x$Preprocess$Original)~rep(1:x$Preprocess$Variates, each=x$Preprocess$Iterations), ...)
	temp <- (1:x$Preprocess$Variates)
	for(i in 1:nrow(x$Observation)) {
		lines(x$Observation[i,]~temp, col="red", type="l")
	}
}
