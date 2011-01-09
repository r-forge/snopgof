gof.optimize <-
function(null, alternate, typeOneError=.05, weights=NULL, verbose=FALSE, iterations=50) {
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
		cat("  >> ",objective,"\n")
		objective # * -1
	}
	cat("Optimizing over weights. Power at alpha = ", typeOneError, ":\n")
	r <- optim(weights,fn,method="SANN", control=list(maxit=iterations))
	matrix(r$par, nrow=null$Variates)
}

