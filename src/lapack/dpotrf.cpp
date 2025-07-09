#include "newscf/lapackinterface.hpp"

extern "C" {
    void dpotrs_(char* uplo, int* n, int* nrhs,
                  double* A, int* lda,
                  double* B, int* ldb, int* info);
}

void newscf_dpotrs(char* uplo, int* n, int* nrhs,
                  double* A, int* lda,
                  double* B, int* ldb, int* info) {
    dpotrs_ (uplo, n, nrhs, A, lda, B, ldb, info);
}