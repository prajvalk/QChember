#include "newscf/lapackinterface.hpp"

extern "C" {
	void dtrsm_ (char* side, char* uplo, char* trans, char* diag,
  int* m, int* n, double* alpha,
  double* A, int* lda, double* B, int* ldb);
}

void newscf_dtrsm (char* side, char* uplo, char* trans, char* diag,
int* m, int* n, double* alpha,
double* A, int* lda, double* B, int* ldb) {
	dtrsm_(side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
}