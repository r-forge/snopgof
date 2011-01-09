gof <-
function(observed, simulated, weights=NULL) {
	if (is.null(weights)) {
		weights = matrix(sapply(1:simulated$Variates, 
			function(i) sapply(1:simulated$Variates, 
			function(j) 
				if (i==j) {
					1/simulated$Variates
				} else {
					0
				} 
			)), ncol=simulated$Variates)
	}
	if (class(observed) != "numeric" | class(simulated) != "gof.preprocess") {
		stop("Parameters are invalid.")
	}
	if (length(observed) != simulated$Variates) {
		stop("Dimensionality of function parameters do not match.")
	}
	tailComparisons.observed = gof.compare(observed, simulated$Original)
	tailComparisons.mean = apply(simulated$Comparisons,1,mean)
	tailComparisons.observed.centered = tailComparisons.observed - tailComparisons.mean
	tailComparisons.simulated.centered = apply(simulated$Comparisons, 2, function(i) i - tailComparisons.mean)
	testStatistic.observed = as.numeric(tailComparisons.observed.centered %*% weights %*% tailComparisons.observed.centered)
	testStatistic.simulated = as.numeric(sapply(1:simulated$Iterations, function(i) tailComparisons.simulated.centered[,i] %*% weights %*% tailComparisons.simulated.centered[,i]))
	#testStatistic.observed = as.numeric(tailComparisons.observed %*% weights %*% tailComparisons.observed)
	#testStatistic.simulated = as.numeric(sapply(1:simulated$Iterations, function(i) simulated$Comparisons[,i] %*% weights %*% simulated$Comparisons[,i]))
	edf.observed = sum(testStatistic.simulated <= testStatistic.observed) / simulated$Iterations
	p <- 1 - abs(1 - 2 * edf.observed)
	return = list(p=p, v.obs=testStatistic.observed, v.sim=testStatistic.simulated)
	class(return) <- "gof.result"
	return
}

