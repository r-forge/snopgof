gof.optimize <-
function(null, alternate, typeOneError=.05, weights=NULL, 
		verbose=FALSE, sannIterations=50, permIterations=127) {
	if (null$Variates < 2) {
		1
	} else {
		if (is.null(weights)) {
			weights = matrix(sapply(1:null$Variates, 
							function(i) sapply(1:null$Variates, 
										function(j) 
											if (i==j) {
												1/null$Variates
											} else {
												0
											} 
								)), ncol=null$Variates)
		}
		if (class(null) != "gof.preprocess" | class(alternate) != "matrix" | typeOneError <= 0 | typeOneError >= 1 | class(weights) != "matrix") {
			stop("Parameters are invalid.")
		}
		if (null$Variates != ncol(alternate) | null$Variates != nrow(weights) | null$Variates != ncol(weights)) {
			stop("Variates do not match in null and alternate")
		}
		fn <- function(w) {
			w <- matrix(w, nrow=null$Variates)
			a <- sapply(1:nrow(alternate), function(i) gof(alternate[i,],null,w)$p)
			p = sum(a <= typeOneError) / nrow(alternate)
			if(verbose) {
				cat("Testing power, summary of p-values for weighting ",w,":\n ")
				cat(" w. alpha .01 = ", sum(a < .01) / nrow(alternate), "\n",
						" w. alpha .05 = ", sum(a < .05) / nrow(alternate), "\n",
						" w. alpha .10 = ", sum(a < .10) / nrow(alternate), "\n",
						" w. alpha .25 = ", sum(a < .25) / nrow(alternate), "\n")
			}
			objective = p
			cat("  >> ",objective," (", diag(w) ,")\n")
			objective * -1
		}
		mostPowerfulWeights <- NULL
		mostPowerfulWeights.p = 0
		cat("Testing permutations of diagonals:\n  >> POWER  ( Diag. of candidate weight )\n")
		# Find the depth that we can go in combinations before exceeding
		# the permutation iterations threshold:
		maxLevel = 0
		count = 0
		while ((count < permIterations) & (maxLevel < null$Variates)) {
			count = count + choose(null$Variates, null$Variates-maxLevel)
			maxLevel = maxLevel+1
		}
		if (count > permIterations) {
			maxLevel = maxLevel-1
			count = count - choose(null$Variates, null$Variates-maxLevel)
		}
		## Create the set of weightings:
		makeMatrixCandidate <- function(vars) {
			vars = vars[vars!=0]
			if ( length(vars) != length(unique(vars)) | length(intersect(vars,1:null$Variates)) != length(vars)) {
				stop("Bad input to make matrix candidate")
			}
			X <- matrix(0, null$Variates, null$Variates)
			if (length(vars) > 0) {
				diag(X)[vars] = 1/length(vars)
			} else {
				diag(X) = 1/null$Variates
			}
			X
		}
		instructions <- matrix(unlist(sapply(0:(maxLevel-1), function(i) {
									(X <- combn(null$Variates,null$Variates-i))
									sapply(1:ncol(X), function(j) {  c(X[,j], rep(0, null$Variates-nrow(X))) } )
								})),nrow=null$Variates)
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

