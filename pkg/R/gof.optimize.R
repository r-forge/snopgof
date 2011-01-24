gof.optimize <-
function(null, alternate, typeOneError=.05, weights=NULL, 
		verbose=FALSE, sannIterations=1000, permIterations=1023) {
	dimensions = null$Variates + 1
	if (null$Variates < 2) {
		1
	} else {
		if (is.null(weights)) {
			weights = matrix(0, dimensions, dimensions)
			diag(weights) <- 1 / (dimensions)
		}
		if (class(null) != "gof.preprocess" | class(alternate) != "matrix" | typeOneError <= 0 | typeOneError >= 1 | class(weights) != "matrix") {
			stop("Parameters are invalid.")
		}
		if (null$Variates != ncol(alternate) | dimensions != nrow(weights) | dimensions != ncol(weights)) {
			stop("Variates do not match in null and alternate")
		}
		fn <- function(w) {
			w <- matrix(w, nrow=dimensions)
			a <- gof(alternate, null, w)$p
			p = sum(a <= typeOneError) / nrow(alternate)
			if(verbose) {
				cat("Testing power, summary of p-values for weighting ",w,":\n ")
				cat(" w. alpha .01 = ", sum(a < .01) / nrow(alternate), "\n",
						" w. alpha .05 = ", sum(a < .05) / nrow(alternate), "\n",
						" w. alpha .10 = ", sum(a < .10) / nrow(alternate), "\n",
						" w. alpha .25 = ", sum(a < .25) / nrow(alternate), "\n")
			}
			objective = p
			if (verbose) {
				cat("  >> ",objective," (", diag(w) ,")\n")
			}
			objective * -1
		}
		mostPowerfulWeights <- NULL
		mostPowerfulWeights.p = 0
		cat("Testing permutations of diagonals:\n  >> POWER  ( Diag. of candidate weight )\n")
		# Find the depth that we can go in combinations before exceeding
		# the permutation iterations threshold:
		maxLevel = 0
		count = 0
		while ((count < permIterations) & (maxLevel < dimensions)) {
			count = count + choose(dimensions, dimensions-maxLevel)
			maxLevel = maxLevel+1
		}
		if (count > permIterations) {
			maxLevel = maxLevel-1
			count = count - choose(dimensions, dimensions-maxLevel)
		}
		## Create the set of weightings:
		makeMatrixCandidate <- function(vars) {
			vars = vars[vars!=0]
			if ( length(vars) != length(unique(vars)) | length(intersect(vars,1:dimensions)) != length(vars)) {
				stop("Bad input to make matrix candidate")
			}
			X <- matrix(0, dimensions, dimensions)
			if (length(vars) > 0) {
				diag(X)[vars] = 1/length(vars)
			} else {
				diag(X) = 1/dimensions
			}
			X
		}
		instructions <- matrix(unlist(sapply(0:(maxLevel-1), function(i) {
									(X <- combn(dimensions,dimensions-i))
									sapply(1:ncol(X), function(j) {  c(X[,j], rep(0, dimensions-nrow(X))) } )
								})),nrow=dimensions)
		candidate.results <- sapply(1:ncol(instructions), function(i) {
					fn(makeMatrixCandidate(instructions[,i]))
				})
		mostPowerfulWeights <- makeMatrixCandidate(instructions[,which(candidate.results == min(candidate.results))[1]])
		mostPowerfulWeights.p <- min(candidate.results)
		
		cat("Optimizing over weights. Power at alpha = ", typeOneError, ":\n")
		r <- optim(mostPowerfulWeights, fn,method="SANN", control=list(maxit=sannIterations))
		r$par
	}
}

