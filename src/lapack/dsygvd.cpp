#include "newscf/lapackinterface.hpp"

extern "C" {
    void dsygvd_(int* itype, char* jobz, char* uplo,
                 int* n, double* A, int* lda,
                 double* B, int* ldb,
                 double* w, double* work, int* lwork,
			     int* iwork, int* liwork,
                 int* info);
}

void newscf_dsygvd(int* itype, char* jobz, char* uplo,
                 int* n, double* A, int* lda,
                 double* B, int* ldb,
                 double* w, double* work, int* lwork,
                 int* iwork, int* liwork,
                 int* info) {
    dsygvd_ (itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, iwork, liwork, info);
}