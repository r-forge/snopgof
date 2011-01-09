gof.preprocess <-
function(simulated) {
	if (class(simulated) != "matrix") {
		stop("Invalid input.")
	}
	simulations=nrow(simulated)
	tailComparisons.simulated = sapply(1:simulations, function(i) gof.compare(simulated[i,], simulated))
	ret <- list(Comparisons=tailComparisons.simulated, Iterations=dim(simulated)[1], Variates=dim(simulated)[2], Original=simulated)
	class(ret) <- "gof.preprocess"
	ret
}

