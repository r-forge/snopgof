gof <-
function(observed, simulated, weights=NULL) {
	if (class(observed) == "numeric") {
		observed <- matrix(observed,nrow=1)
	}

	if (class(observed) != "matrix" | class(simulated) != "gof.preprocess") {
		stop("simulated must be a gof.preprocess object.")
	}
	if (ncol(observed) != simulated$Variates) {
		stop("Dimensionality of function parameters do not match.")
	}
	observations = nrow(observed)
	
	if (is.null(weights)) {
		if (length(observed) > 1) {
			weights = matrix(0, simulated$Variates+1, simulated$Variates+1)
			diag(weights) <- 1 / (simulated$Variates+1)
		} else {
			weights = 1
		}
	}
	
	foldedObservations <- sapply(1:observations, function(i) abs(observed[i,]-simulated$Expectation))
	p <- .Call("doTest", as.matrix(foldedObservations), as.matrix(simulated$FoldedSimulations),
				as.matrix(simulated$CenteredComparisons), as.matrix(simulated$ComparisonCrease), 
								as.matrix(weights))
	
	#	Comparisons in R instead of C
	#	tailComparisons.observed = sapply(1:observations, function (i) gof.compare(abs(observed[i,]-simulated$Expectation), simulated$FoldedSimulations) )
	#	tailComparisons.observed.centered = abs(tailComparisons.observed - simulated$ComparisonCrease)
	#	testStatistic.observed =  as.numeric(sapply(1:observations, function(i) t(tailComparisons.observed.centered[,i]) %*% weights %*% tailComparisons.observed.centered[,i]))
	#	testStatistic.simulated = as.numeric(sapply(1:simulated$Iterations, function(i) t(simulated$CenteredComparisons[,i]) %*% weights %*% simulated$CenteredComparisons[,i]))
	#	edf.observed = sapply(1:observations, function (i) sum(testStatistic.simulated <= testStatistic.observed[i]) / simulated$Iterations )
	#	p <- 1 - abs(1 - 2 * edf.observed)
	mhd <- sapply(1:observations, function(i) t(observed[i,] - simulated$Expectation) %*% simulated$InverseCovariance %*% (observed[i,] - simulated$Expectation) )
	p.mhd <- sapply(1:observations, function (i) 1 - abs(1 - 2 * sum(mhd[i] >= simulated$MHDistances)/simulated$Iterations) )
	return = list(p=p, #v.obs=testStatistic.observed, v.sim=testStatistic.simulated,
			Observation=observed, Preprocess=simulated, MHDistance=mhd, p.MHD = p.mhd)
	class(return) <- "gof.result"
	return
}

print.gof.result <- function (x, ...) {
	with(x, cat("snopgof results\nVariates:", Preprocess$Variates ,"\nSimulations:", Preprocess$Iterations ,"\nP-Value:", p,"\nP-Value(MHD):", p.MHD))
}

plot.gof.result <- function (x, standardize=3, violin=TRUE, 
		main="SNOPGOF Goodness of Fit Test", ylim=NULL, xlim=NULL, 
		xlab=NULL, ylab=NULL, key=NULL, perc=.05, ...) {
	itns <- x$Preprocess$Iterations
	vars <- x$Preprocess$Variates
	sims <- x$Preprocess$Original
	obs <- x$Observation
	n.obs <- nrow(obs)
	
	if (standardize==3) {
		sims.median <- apply(sims, 2, median)
		sims.min <- apply(sims, 2, min)
		sims.max <- apply(sims, 2, max)
		sims <- sapply(1:ncol(sims), function(i) (sims[,i] - sims.median[i])/(sims.max[i] - sims.min[i] ) )
		obs <- matrix(sapply(1:ncol(sims), function(i) (obs[,i] - sims.median[i])/(sims.max[i] - sims.min[i] ) ), nrow=n.obs )
	} else if (standardize==2) {
		sims.min <- apply(sims, 2, min)
		sims.max <- apply(sims, 2, max)
		sims <- sapply(1:ncol(sims), function(i) (sims[,i] - sims.min[i])/(sims.max[i] - sims.min[i] ) )
		obs <- matrix(sapply(1:ncol(sims), function(i) (obs[,i] - sims.min[i])/(sims.max[i] - sims.min[i] ) ), nrow=n.obs )
	} else if (standardize==1) {
		sims.mean <- apply(sims, 2, mean)
		sims <- sapply(1:ncol(sims), function(i) (sims[,i] - sims.mean[i]) )
		obs <- matrix(sapply(1:ncol(sims), function(i) (obs[,i] - sims.mean[i]) ), nrow=n.obs )
	}
	
	if (is.null(ylim)) {
		ylim = c(min(obs, sims), max(obs, sims))
	}
	if (is.null(xlim)) {
		xlim = c(0, ncol(obs)+1)
	}
	if (is.null(xlab)) {
		xlab= paste( paste("p:", round(x$p, 3), collapse = " "), paste("p(MHD):", round(x$p.MHD, 3), collapse = " "), collapse = "\n")
	}
	if (is.null(ylab)) {
		ylab = "Statistic Values"
	}
	xAxis <- (1:vars)

	plot(obs[1,]~xAxis, col="white", type="p",  ylim=ylim, xlim=xlim, main=main, xlab=xlab, ylab=ylab, axes=FALSE, ...)
	if (!is.null(key)) {
		if (length(key) != ncol(obs)) {
			stop("Key length does not match the number of variates.")
		}
		axis(1, at=xAxis, lab=key)
	} else {
		axis(1, at=xAxis, lab=paste("v", xAxis, sep=""))
	}

	ind.lower = round(itns * perc/2)
	ind.upper = round(itns * (1-perc/2))
	yperc.lower = sapply(1:ncol(sims), function(i)  sort(sims[,i])[ind.lower]  )
	yperc.upper = sapply(1:ncol(sims), function(i)  sort(sims[,i])[ind.upper]  )
	lines(yperc.lower~xAxis, lty=3, col = "gray", lwd=3)
	lines(yperc.upper~xAxis, lty=3, col = "gray", lwd=3)

	if (violin) {
		require(vioplot)
		for (i in 1:ncol(sims)) {
			vioplot(sims[,i], at=xAxis[i], add=TRUE, col="gray", wex=.75, ...)
		}
	} else {
		boxplot(as.numeric(sim)~rep(1:vars, each=itns), add=TRUE, ...)
	}
	for(i in 1:nrow(obs)) {
		lines(obs[i,]~xAxis, col="red", type="l", lwd=1, ...)
		lines(obs[i,]~xAxis, col="red", type="p", lwd=3, pch=19, ...)
		text(xAxis, obs[i,], labels=round(x$Observation[i,],3), pos=4)
	}

}
