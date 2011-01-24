#include <R.h>
#include <Rdefines.h>
#include <math.h>

SEXP comparePreprocess(SEXP simulations) {
	SEXP Rdim = getAttrib(simulations, R_DimSymbol);
	SEXP result;
	int nSimulations = INTEGER(Rdim)[0];
	int nVariates = INTEGER(Rdim)[1];
	PROTECT(result = allocMatrix(REALSXP, nVariates+1, nSimulations));
	for (int simIndex1=0; simIndex1<nSimulations; simIndex1++) {
		for (int threshold=0; threshold<=nVariates; threshold++) {
			REAL(result)[simIndex1 + threshold*nSimulations]=0;
		}
	}
	for (int simIndex1=0; simIndex1<nSimulations; simIndex1++) {
		for (int simIndex2=simIndex1+1; simIndex2<nSimulations; simIndex2++) {
			int num1geq2 =0, num2geq1=0;
			for (int varIndex=0; varIndex<nVariates; varIndex++) {
				if ( REAL(simulations)[simIndex1 + varIndex*nSimulations] >
					REAL(simulations)[simIndex2 + varIndex*nSimulations] ) {
					num1geq2++;
				} else if ( REAL(simulations)[simIndex1 + varIndex*nSimulations] <
					REAL(simulations)[simIndex2 + varIndex*nSimulations] ) {
					num2geq1++;
				} else {
					num1geq2++;
					num2geq1++;
				}
			}
			REAL(result)[simIndex1*(nVariates+1) + num1geq2]++;
			REAL(result)[simIndex2*(nVariates+1) + num2geq1]++;
		}
	}
	UNPROTECT(1);
	return result;
}
// .Call("doTest", foldedObservations, simulated$FoldedSimulations, simulated$FoldedComparisons, simulated$ComparisonCrease, weights)
SEXP doTest(SEXP observations, SEXP simulations, SEXP comparisons, SEXP comparisonCrease, SEXP weights) {
	SEXP Rdim;
	Rdim = getAttrib(observations, R_DimSymbol);
	int nVariates = INTEGER(Rdim)[0];
	int nObservations = INTEGER(Rdim)[1];
	Rdim = getAttrib(simulations, R_DimSymbol);
	int nSimulations = INTEGER(Rdim)[0];
	SEXP result;
	PROTECT(result = allocVector(REALSXP, nObservations));
//	Rprintf("Obs: %u, Var: %u, Sim: %u\n", nObservations, nVariates, nSimulations);

	// Find the null distribution of the omnibus test statistic:
	double *omnibusStatisticNull = (double*)malloc(sizeof(double) * nSimulations);
	for (int simIndex=0; simIndex<nSimulations; simIndex++) {
		omnibusStatisticNull[simIndex] = 0;
		for (int i=0; i<=nVariates; i++) {
			for (int j=0; j<=nVariates; j++) {
				omnibusStatisticNull[simIndex] +=
						// Not sure which order to put these in!
//						REAL(comparisons)[i*nSimulations + simIndex] *
//						REAL(comparisons)[j*nSimulations + simIndex] *
						REAL(comparisons)[i + simIndex * (nVariates+1)] *
						REAL(comparisons)[j + simIndex * (nVariates+1)] *
						REAL(weights)[i + (nVariates+1)*j];
			}
		}
//		Rprintf("omnibusStat %u: %f\n", simIndex, omnibusStatisticNull[simIndex]);
	}

	// Initialize a vector for storing tail statistics for each observation
	double *tailStatistic = (double*)malloc(sizeof(double)*(nVariates+1));

	for (int obsIndex=0; obsIndex<nObservations; obsIndex++) {
		// Reset the tail statistics to minus the crease (saves time later)
		for (int i=0; i<=nVariates; i++) {
			tailStatistic[i] = -1 * REAL(comparisonCrease)[i];
		}
		// Do tail comparisons between this observation and the simulated data.
		for (int simIndex=0; simIndex<nSimulations; simIndex++) {
			int obsGeqSim = 0;
			// How many coordinates is the observation >= simulation for?
			for (int varIndex=0; varIndex<nVariates; varIndex++) {
//				if ( REAL(observations)[obsIndex + varIndex * nObservations] >= REAL(simulations)[simIndex + varIndex * nSimulations] ) {
				if ( REAL(observations)[obsIndex * nVariates + varIndex] >= REAL(simulations)[simIndex + varIndex * nSimulations] ) {
					obsGeqSim++;
				}
			}
			// Increment that tail statistic.
			tailStatistic[obsGeqSim]++;
		}
		// Take the absolute value, completing the folding
		for (int i=0; i<=nVariates; i++) {
			tailStatistic[i] = tailStatistic[i];
		}
		// Now we have calculated the tail statistic for the observation. Calculate the omnibus test statistic:
		double omnibus = 0;
		for (int i=0; i<=nVariates; i++) {
			for (int j=0; j<=nVariates; j++) {
				omnibus += tailStatistic[i] * tailStatistic[j] * REAL(weights)[i + (nVariates+1)*j];
			}
		}
		// Find the p-value:
		REAL(result)[obsIndex]=0;
		for (int simIndex=0; simIndex < nSimulations; simIndex++) {
			if (omnibusStatisticNull[simIndex] <= omnibus) {
				REAL(result)[obsIndex]++;
			}
		}
		REAL(result)[obsIndex] = 1 - fabs( 1 - 2 * REAL(result)[obsIndex] / nSimulations );
	}
	// Free the vector of tail statistics
	free(tailStatistic);
	free(omnibusStatisticNull);
	UNPROTECT(1);
	return result;
}
