#include "newscf/lapackinterface.hpp"

extern "C" {
    void dorgqr_(int* m, int* n,
                  double* A, int* lda,
                  double* tau, double* work, int* lwork,
                  int* info);
}

void newscf_dorgqr(int* m, int* n, int* k,
				  double* A, int* lda,
				  double* tau, double* work, int* lwork,
				  int* info) {
    dorgqr_ (m, n, A, lda, tau, work, lwork, info);
}