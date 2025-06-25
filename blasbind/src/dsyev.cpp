#include "lapackinterface.hpp"

extern "C" {
    void dsyev_(char* jobz, char* uplo, int* n,
                 double* A, int* lda, double* w,
                 double* work, int* lwork, int* info);
}

void newscf_dsyev(char* jobz, char* uplo, int* n,
                double* A, int* lda, double* w, double* work,
                int* lwork, int* info) {
    dsyev_ (jobz, uplo, n, A, lda, w, work, lwork, info);
}