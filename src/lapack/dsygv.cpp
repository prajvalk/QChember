#include "newscf/lapackinterface.hpp"

extern "C" {
    void dsygv_(int* itype, char* jobz, char* uplo,
                 int* n, double* A, int* lda,
                 double* B, int* ldb,
                 double* w, double* work, int* lwork,
                 int* info);
}

void newscf_dsygv(int* itype, char* jobz, char* uplo,
                 int* n, double* A, int* lda,
                 double* B, int* ldb,
                 double* w, double* work, int* lwork,
                 int* info) {
    dsygv_ (itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, info);
}