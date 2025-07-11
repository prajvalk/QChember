#ifndef NEWSCF_LAPACKINTERFACE_HPP
#define NEWSCF_LAPACKINTERFACE_HPP

// --- BLAS routines ---
void newscf_dgemm (char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* alpha,
    double* A, int* LDA,double* B, int* LDB, double* beta, double* C, int* LDC);

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

void newscf_dsygvd(int* itype, char* jobz, char* uplo,
                 int* n, double* A, int* lda,
                 double* B, int* ldb,
                 double* w, double* work, int* lwork,
                 int* iwork, int* liwork,
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

// Symmetric Linsolve Routines

void newscf_dsysv (char* uplo, int* n, int* nrhs, double* a, int* lda,
                   int* ipiv, double* b, int* ldb, double* work, int* lwork, int* info);

#include "common.hpp"
using newscf::handle_error;

// TODO: Instead of enum style backends do some C++20 enum-templates
// TODO: for highly optimized structures
// current logic requires all backends to be available for link at compile time
enum GeneralGEV_Backend {
  DSYGV,  // QR
  DSYGVD, // Divide and Conquer
};

// General GEV Eigensolver Wrapper
class GeneralGEV {
public:
  GeneralGEV_Backend backend;

  // Used by DSYGxx
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

  // Used by only DSYGD
  int* IWORK;
  int  LIWORK;

  GeneralGEV() {
    backend = DSYGVD;
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

    IWORK = nullptr;
    LIWORK = -1;
  }

  ~GeneralGEV() {
    delete[] W;
    delete[] A;
    delete[] B;
    delete[] WORK;
    delete[] IWORK;
  }

  inline void set_itype (const int type) {
    ITYPE = type;
  }

  inline void compute_eigenvectors (bool eigvec) {
    JOBZ = eigvec ? 'V' : 'U';
  }

  inline void set_backend (const GeneralGEV_Backend b) {
    backend = b;
  }

  inline void set_N (const int size) {
    N = size;
    LDA = size;
    LDB = size;
  }

  inline int init() {
    LWORK = -1;
    WORK = new double[1];
    if (backend == DSYGVD) {
      LIWORK = -1;
      IWORK = new int[1];
    }
    if (backend == DSYGV) newscf_dsygv (&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, &INFO);
    else if (backend == DSYGVD) newscf_dsygvd (&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, IWORK, &LIWORK, &INFO);

    if (INFO == 0) {
      LWORK = static_cast<int>(WORK[0]);
      delete[] WORK;
      WORK = new double[LWORK];
      A = new double[N * N];
      B = new double[N * N];
      W = new double[N];

      if (backend == DSYGVD) {
        LIWORK = IWORK[0];
        delete[] IWORK;
        IWORK = new int[LIWORK];
      }
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
    if (backend == DSYGV) newscf_dsygv (&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, &INFO);
    else if (backend == DSYGVD) newscf_dsygvd (&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, IWORK, &LIWORK, &INFO);
    if (INFO != 0)
      HANDLE_ERROR("GeneralGEV (backend="+std::to_string(backend)+") failed with INFO = "+ std::to_string(INFO), 90);
    return INFO;
  }

  inline void copy_eigenvectors (double* result) {
    newscf::memcpy(result, A, LDA * N * sizeof(double));
  }
};

enum GeneralLinsolve_Backend {
  DSYSV // LDL + Rook v1/Bunch–Kaufman (depends on LAPACK version)
  // TODO: DSYSV_ROOK // LDL + Rook v2 Pivoting
  // TODO: DSYSV_RK // LDL + Rook v3 Pivoting
  // TODO: DSYSV_AA // Aasen
  // TODO: DSYSV_AA_2 // Aasen 2 Stage
};

class GeneralLinsolve {
public:
  GeneralLinsolve_Backend backend;

  // used by all DSYSVxx
  char  UPLO;
  int   N;
  int   NRHS;
  double* A;
  int   LDA;
  int* IPIV;
  double* B;
  int   LDB;
  double* WORK;
  int   LWORK;
  int INFO;

  GeneralLinsolve() {
    backend = DSYSV;
    UPLO = 'U';
    N = 1;
    NRHS = 1;
    LDA = 1;
    LDB = 1;
    LWORK = -1;
    INFO = 0;
    WORK = nullptr;
    A = nullptr;
    B = nullptr;
    IPIV = nullptr;
  }

  ~GeneralLinsolve() {
    delete[] WORK;
    delete[] A;
    delete[] B;
    delete[] IPIV;
  }

  void set_uplo (const char uplo) {
    UPLO = uplo;
  }

  void set_N (const int n) {
    N = n;
    LDA = n;
    LDB = n;
  }

  void set_NRHS (const int nrhs) {
    NRHS = nrhs;
  }

  void set_backend (const GeneralLinsolve_Backend b) {
    backend = b;
  }

  int init() {
    LWORK = -1;
    INFO = 0;
    WORK = new double[1];
    if (backend == DSYSV) {
      newscf_dsysv (&UPLO, &N, &NRHS, A, &LDA, IPIV, B, &LDB, WORK, &LWORK, &INFO);
    }
    if (INFO == 0) {
      LWORK = static_cast<int>(WORK[0]);
      delete[] WORK;
      WORK = new double[LWORK];
      A = new double[N * N];
      B = new double[N * NRHS];
      IPIV = new int[N];
    }
    return INFO;
  }

  void set_A (const double* rA) {
    newscf::memcpy(A, rA, LDA * N * sizeof(double));
  }

  void set_B (const double* rB) {
    newscf::memcpy(B, rB, LDB * NRHS * sizeof(double));
  }

  int run() {
    if (backend == DSYSV) newscf_dsysv (&UPLO, &N, &NRHS, A, &LDA, IPIV, B, &LDB, WORK, &LWORK, &INFO);
    if (INFO != 0) {
      HANDLE_ERROR("GeneralLinsolve (backend="+std::to_string(backend)+") failed with INFO = "+ std::to_string(INFO), 90);
    }
    return INFO;
  }

  void copy_result (double* result) {
    newscf::memcpy(result, B, LDB * NRHS * sizeof(double));
  }
};

#endif