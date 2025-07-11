#include "newscf/lapackinterface.hpp"

extern "C" {
	void dsysv_ (char* uplo, int* n, int* nrhs, double* a, int* lda,
					   int* ipiv, double* b, int* ldb, double* work, int* lwork, int* info);
}

void newscf_dsysv (char* uplo, int* n, int* nrhs, double* a, int* lda,
				   int* ipiv, double* b, int* ldb, double* work, int* lwork, int* info) {
	dsysv_ (uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}