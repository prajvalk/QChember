#include "newscf/lapackinterface.hpp"

extern "C" {
    void dsyevd_(char* jobz, char* uplo, int* n,
                  double* A, int* lda, double* w,
                  double* work, int* lwork,
                  int* iwork, int* liwork,
                  int* info);
}

void newscf_dsyevd(char* jobz, char* uplo, int* n,
                  double* A, int* lda, double* w,
                  double* work, int* lwork,
                  int* iwork, int* liwork,
                  int* info) {
    dsyevd_ (jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info);
}