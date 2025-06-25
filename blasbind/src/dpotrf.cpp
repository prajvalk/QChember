#include "lapackinterface.hpp"

extern "C" {
    void _dpotrs(char* uplo, int* n, int* nrhs,
                  double* A, int* lda,
                  double* B, int* ldb, int* info);
}

void newscf_dpotrs(char* uplo, int* n, int* nrhs,
                  double* A, int* lda,
                  double* B, int* ldb, int* info) {
    _dpotrs (uplo, n, nrhs, A, lda, B, ldb, info);
}