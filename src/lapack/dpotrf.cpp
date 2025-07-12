#include "newscf/lapackinterface.hpp"

extern "C" {
	void dpotrf_ (char* uplo, int* n,
				  double* A, int* lda, int* info);
}

void newscf_dpotrf(char* uplo, int* n,
				  double* A, int* lda, int* info) {
	dpotrf_(uplo, n, A, lda, info);
}