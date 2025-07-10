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

#include "common.hpp"
using newscf::handle_error;

// DSYGV
class DSYGV {
public:
  int     ITYPE;
  char    JOBZ;
  char    UPLO;
  int     N;
  int     LDA;
  int     LDB;
  int     LWORK;
  int     INFO;
  double* W;
  double* WORK;
  double* A;
  double* B;

  DSYGV() {
    ITYPE = 1;
    JOBZ = 'V';
    UPLO = 'U';
    N = 1;
    LDA = 1;
    LDB = 1;
    LWORK = -1;
    INFO = 0;
    W = nullptr;
    WORK = nullptr;
    A = nullptr;
    B = nullptr;
  }

  ~DSYGV() {
    delete[] W;
    delete[] A;
    delete[] B;
    delete[] WORK;
  }

  inline void set_itype (const int type) {
    ITYPE = type;
  }

  inline void compute_eigenvectors (bool eigvec) {
    JOBZ = eigvec ? 'V' : 'U';
  }

  inline void set_N (const int size) {
    N = size;
    LDA = size;
    LDB = size;
  }

  inline int init() {
    LWORK = -1;
    WORK = new double[1];
    newscf_dsygv (&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, &INFO);
    if (INFO == 0) {
      LWORK = static_cast<int>(WORK[0]);
      delete[] WORK;
      WORK = new double[LWORK];
      A = new double[N * N];
      B = new double[N * N];
      W = new double[N];
    }
    return INFO;
  }

  inline void set_A (const double* rA) {
    newscf::memcpy(A, rA, LDA * N * sizeof(double));
  }

  inline void set_B (const double* rB) {
    newscf::memcpy(B, rB, LDB * N * sizeof(double));
  }

  inline int run() {
    newscf_dsygv (&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, &INFO);
    if (INFO != 0)
      HANDLE_ERROR("DSYGV failed with INFO = "+ std::to_string(INFO), 90);
    return INFO;
  }

  inline void copy_eigenvectors (double* result) {
    newscf::memcpy(result, A, LDA * N * sizeof(double));
  }
};



#endif