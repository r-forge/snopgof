gof.mhd.asymptotic <-
function(observed, null.mean, null.cov, twoTailed = FALSE) {
	if (class(observed) != "numeric" | class(null.mean) != "numeric" | class(null.cov) != "matrix") {
		stop("Parameters are invalid.")
	}
	if (length(observed) != length(null.mean) | length(observed) != nrow(null.cov) | length(observed) != ncol(null.cov)) {
		stop("Dimensionality of function parameters do not match.")
	}
	tstat <- (observed - null.mean) %*% solve(null.cov) %*% (observed - null.mean)
	if (twoTailed) {
		pval <- 1- abs(1 - 2 * pchisq(tstat, length(observed)))
	} else {
		pval <- 1 - pchisq(tstat, length(observed))
	}
	
	return = list(p=pval, v.obs=tstat)
	class(return) <- "gof.mhd.asymptotic.result"
	return
}
