#include "lapackinterface.hpp"

extern "C" {
    void dgesv_(int* n, int* nrhs,
                 double* A, int* lda,
                 int* ipiv, double* B, int* ldb,
                 int* info);
}

void newscf_dgesv(int* n, int* nrhs,
                 double* A, int* lda,
                 int* ipiv, double* B, int* ldb,
                 int* info) {
    dgesv_ (n, nrhs, A, lda, ipiv, B, ldb, info);
}