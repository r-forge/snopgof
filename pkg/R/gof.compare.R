gof.compare <-
function(observed, simulated) {
	if (class(observed) != "numeric" | class(simulated) != "matrix") {
		stop("Parameters are invalid.")
	}
	if (length(observed) != dim(simulated)[2]) {
		stop("Dimensionality of function parameters do not match.")
	}
	variates = length(observed)
	simulations = dim(simulated)[1]
	a <- sapply(0:variates, function (threshold) {sum(sapply(1:simulations, 
		function (i) sum(observed >= simulated[i,]) == threshold)) 
	})
	class(a) <- "gof.comparison"
	a
}
