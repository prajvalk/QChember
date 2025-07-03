#ifndef NEWSCF_LAPACKINTERFACE_HPP
#define NEWSCF_LAPACKINTERFACE_HPP

/*
// --- Level 3 BLAS ---
void newscf_dgemm_(char* transa, char* transb,
                 int* m, int* n, int* k,
                 double* alpha, double* A, int* lda,
                 double* B, int* ldb,
                 double* beta, double* C, int* ldc);

void newscf_sgemm_(char* transa, char* transb,
                 int* m, int* n, int* k,
                 float* alpha, float* A, int* lda,
                 float* B, int* ldb,
                 float* beta, float* C, int* ldc);

// --- Level 2 BLAS ---
void newscf_dgemv_(char* trans, int* m, int* n,
                 double* alpha, double* A, int* lda,
                 double* x, int* incx,
                 double* beta, double* y, int* incy);

void newscf_sgemv_(char* trans, int* m, int* n,
                 float* alpha, float* A, int* lda,
                 float* x, int* incx,
                 float* beta, float* y, int* incy);

// --- Level 1 BLAS ---
void newscf_daxpy_(int* n, double* alpha,
                 double* x, int* incx,
                 double* y, int* incy);

void newscf_saxpy_(int* n, float* alpha,
                 float* x, int* incx,
                 float* y, int* incy);

void newscf_dcopy_(int* n, double* x, int* incx,
                 double* y, int* incy);

void newscf_scopy_(int* n, float* x, int* incx,
                 float* y, int* incy);

double newscf_ddot_(int* n, double* x, int* incx, double* y, int* incy);
float  newscf_sdot_(int* n, float* x, int* incx, float* y, int* incy);

void newscf_dscal_(int* n, double* alpha, double* x, int* incx);
void newscf_sscal_(int* n, float* alpha, float* x, int* incx);*/

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