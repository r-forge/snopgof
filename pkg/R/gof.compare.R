gof.compare <-
function(observed, simulated, median=FALSE) {
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
	a <- sapply(0:(variates-1), function (threshold) {sum(sapply(1:simulations, 
		function (i) sum(abs(observed - centroid) > abs(simulated[i,] - centroid)) == threshold)) 
	})
	class(a) <- "gof.comparison"
	a
}

