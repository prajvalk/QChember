#include "lapackinterface.hpp"

extern "C" {
    void _dgeev(char* jobvl, char* jobvr, int* n,
                 double* A, int* lda,
                 double* wr, double* wi,
                 double* vl, int* ldvl,
                 double* vr, int* ldvr,
                 double* work, int* lwork,
                 int* info);
}

void newscf_dgeev(char* jobvl, char* jobvr, int* n,
                 double* A, int* lda,
                 double* wr, double* wi,
                 double* vl, int* ldvl,
                 double* vr, int* ldvr,
                 double* work, int* lwork,
                 int* info) {
    _dgeev (jobvl, jobvr, n, A, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
}