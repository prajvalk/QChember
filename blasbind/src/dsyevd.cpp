#include "lapackinterface.hpp"

extern "C" {
    void _dsyevd(char* jobz, char* uplo, int* n,
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
    _dsyevd (jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info);
}