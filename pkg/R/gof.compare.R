gof.compare <-
function(observed, simulated, centralityMeasure="mean") {
	if (class(observed) != "numeric" | class(simulated) != "matrix") {
		stop("Parameters are invalid.")
	}
	if (length(observed) != dim(simulated)[2]) {
		stop("Dimensionality of function parameters do not match.")
	}
	variates = length(observed)
	simulations = dim(simulated)[1]
	if (centralityMeasure=="median") {
		centralityMeasure = apply(simulated,2,median)
	} else if (centralityMeasure=="mode") {
		centralityMeasure = apply(simulated,2,mode)
	} else if (centralityMeasure=="mean") {
		centralityMeasure = apply(simulated,2,mean)
	} else {
		stop("Invalid parameter for centroid.")
	}
	a <- sapply(0:(variates-1), function (threshold) {sum(sapply(1:simulations, 
		function (i) sum(abs(observed - centralityMeasure) > abs(simulated[i,] - centralityMeasure)) == threshold)) 
	})
	class(a) <- "gof.comparison"
	a
}

