#include "lapackinterface.hpp"

extern "C" {
    void dgeqrf_(int* m, int* n,
                  double* A, int* lda,
                  double* tau, double* work, int* lwork,
                  int* info);
}

void newscf_dgeqrf(int* m, int* n,
                  double* A, int* lda,
                  double* tau, double* work, int* lwork,
                  int* info) {
    dgeqrf_ (m, n, A, lda, tau, work, lwork, info);
}