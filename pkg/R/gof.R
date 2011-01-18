gof <-
function(observed, simulated, weights=NULL) {
	if (is.null(weights)) {
		if (length(observed) > 1) {
			weights = matrix(sapply(1:simulated$Variates, 
							function(i) sapply(1:simulated$Variates, 
										function(j) 
											if (i==j) {
												1/simulated$Variates
											} else {
												0
											} 
								)), ncol=simulated$Variates)
		} else {
			weights = 1
		}
	}
	if (class(observed) != "numeric" | class(simulated) != "gof.preprocess") {
		stop("Parameters are invalid.")
	}
	if (length(observed) != simulated$Variates) {
		stop("Dimensionality of function parameters do not match.")
	}
	tailComparisons.observed = gof.compare(observed, simulated$Original)
	if (class(simulated$Comparisons) == "integer") {
		tailComparisons.mean = mean(simulated$Comparisons)
		tailComparisons.observed.centered = tailComparisons.observed - tailComparisons.mean
		tailComparisons.simulated.centered = simulated$Comparisons - tailComparisons.mean
	} else {
		tailComparisons.mean = apply(simulated$Comparisons,1,mean)
		tailComparisons.observed.centered = tailComparisons.observed - tailComparisons.mean
		tailComparisons.simulated.centered = apply(simulated$Comparisons, 2, function(i) i - tailComparisons.mean)
	}
	testStatistic.observed = as.numeric(tailComparisons.observed.centered %*% weights %*% tailComparisons.observed.centered)
	if(class(tailComparisons.simulated.centered) == "numeric") {
		testStatistic.simulated = as.numeric(sapply(1:simulated$Iterations, function(i) tailComparisons.simulated.centered[i] %*% weights %*% tailComparisons.simulated.centered[i]))
	} else {
		testStatistic.simulated = as.numeric(sapply(1:simulated$Iterations, function(i) tailComparisons.simulated.centered[,i] %*% weights %*% tailComparisons.simulated.centered[,i]))
	}
	edf.observed = sum(testStatistic.simulated <= testStatistic.observed) / simulated$Iterations
	p <- 1 - abs(1 - 2 * edf.observed)
	return = list(p=p, v.obs=testStatistic.observed, v.sim=testStatistic.simulated)
	class(return) <- "gof.result"
	return
}

