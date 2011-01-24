gof.preprocess <-
function(simulated, creasingFunction=median) {
	if (class(simulated) != "matrix") {
		stop("Invalid input.")
	}
	a <- cov(simulated)
	ainv <- solve(a)
	simulations=nrow(simulated)
	expectation = apply(simulated, 2, creasingFunction);
	centeredSimulations <- t(sapply(1:simulations, function(i) simulated[i,] - expectation))
	foldedSimulations = abs(centeredSimulations)
	ttc <- system.time((tailComparisons.simulated <- .Call("comparePreprocess", foldedSimulations)))
	comparisonCrease = apply(tailComparisons.simulated, 1, creasingFunction)
	tailComparisons.simulated.centered = sapply(1:simulations, function (i) tailComparisons.simulated[,i] - comparisonCrease)
	# Only use if C is broken:
	#	tailComparisons.simulated = sapply(1:simulations,
	#			function(i) gof.compare(foldedSimulations[i,], foldedSimulations[-i,]))

	mahalanobisDistances = sapply(1:simulations, function(i) centeredSimulations[i,] %*% ainv %*% centeredSimulations[i,])
	ret <- list(Comparisons=tailComparisons.simulated, CenteredComparisons=tailComparisons.simulated.centered,
			ComparisonCrease = comparisonCrease,
			Iterations=dim(simulated)[1], 
			Variates=dim(simulated)[2], Original=simulated, Covariance=a, InverseCovariance=ainv,
			Expectation=expectation, FoldedSimulations=foldedSimulations,
			CenteredSimulations=centeredSimulations, MHDistances = sort(mahalanobisDistances),
			CreasingFunction=median, ComputeTime=ttc)
	class(ret) <- "gof.preprocess"
	ret
}

print.gof.preprocess <- function (x, ...) {
	with(x, cat("snopgof preprocess object.\nVariates: ", Variates, "\nIterations: ", Iterations, "\nCompute time: ", ComputeTime[[1]]), "s")
}
