#ifndef NEWSCF_LAPACKINTERFACE_HPP
#define NEWSCF_LAPACKINTERFACE_HPP

// --- LAPACK routines ---
void newscf_dsyev(char* jobz, char* uplo, int* n,
                 double* A, int* lda, double* w,
                 double* work, int* lwork, int* info);

void newscf_dgeev(char* jobvl, char* jobvr, int* n,
                 double* A, int* lda,
                 double* wr, double* wi,
                 double* vl, int* ldvl,
                 double* vr, int* ldvr,
                 double* work, int* lwork,
                 int* info);

void newscf_dsyevd(char* jobz, char* uplo, int* n,
                  double* A, int* lda, double* w,
                  double* work, int* lwork,
                  int* iwork, int* liwork,
                  int* info);

void newscf_dgesv(int* n, int* nrhs,
                 double* A, int* lda,
                 int* ipiv, double* B, int* ldb,
                 int* info);

void newscf_dpotrf(char* uplo, int* n,
                  double* A, int* lda, int* info);

void newscf_dpotrs(char* uplo, int* n, int* nrhs,
                  double* A, int* lda,
                  double* B, int* ldb, int* info);

void newscf_dgeqrf(int* m, int* n,
                  double* A, int* lda,
                  double* tau, double* work, int* lwork,
                  int* info);

void newscf_dorgqr(int* m, int* n, int* k,
                  double* A, int* lda,
                  double* tau, double* work, int* lwork,
                  int* info);

void newscf_dsygv(int* itype, char* jobz, char* uplo,
                 int* n, double* A, int* lda,
                 double* B, int* ldb,
                 double* w, double* work, int* lwork,
                 int* info);

void newscf_dggev(char* jobvl, char* jobvr,
                int* n, double* a, int* lda,
                double* b, int* ldb,
                double* alphar, double* alphai,
                double* beta,
                double* vl, int* ldvl,
                double* vr, int* ldvr,
                double* work, int* lwork,
                int* info);



#endif