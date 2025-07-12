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

void newscf_dtrsm (char* side, char* uplo, char* trans, char* diag,
  int* m, int* n, double* alpha,
  double* A, int* lda, double* B, int* ldb);

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
  DSYGVD, // Divide and
  STDTRANS_SYEV
};

class GeneralGEV {
public:
  GeneralGEV_Backend backend;

  int ITYPE;
  char JOBZ;
  char UPLO;
  int N;
  int LDA;
  int LDB;

  int LWORK;  // Total workspace size for STDTRANS_SYEV
  int INFO;

  double* W;      // eigenvalues
  double* WORK;   // combined workspace for STDTRANS_SYEV or DSYGV/DSYGVD
  int* IWORK;     // for DSYGVD
  int LIWORK;

  double* A;
  double* B;

  GeneralGEV()
      : backend(DSYGVD), ITYPE(1), JOBZ('V'), UPLO('U'),
        N(1), LDA(1), LDB(1),
        LWORK(-1), INFO(0),
        W(nullptr), WORK(nullptr), IWORK(nullptr), LIWORK(-1),
        A(nullptr), B(nullptr) {}

  ~GeneralGEV() {
    delete[] W;
    delete[] WORK;
    delete[] IWORK;
    delete[] A;
    delete[] B;
  }

  void set_itype(int type) { ITYPE = type; }
  void compute_eigenvectors(bool eigvec) { JOBZ = eigvec ? 'V' : 'N'; }
  void set_backend(GeneralGEV_Backend b) { backend = b; }
  void set_N(int size) {
    N = size;
    LDA = size;
    LDB = size;
  }

  int init() {
    INFO = 0;

    // Free old buffers if any
    delete[] W; W = nullptr;
    delete[] WORK; WORK = nullptr;
    delete[] IWORK; IWORK = nullptr;
    delete[] A; A = nullptr;
    delete[] B; B = nullptr;

    W = new double[N];
    A = new double[N * N];
    B = new double[N * N];

    if (backend == DSYGVD) {
      LWORK = -1;
      LIWORK = -1;
      WORK = new double[1];
      IWORK = new int[1];

      newscf_dsygvd(&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB,
                    W, WORK, &LWORK, IWORK, &LIWORK, &INFO);
      if (INFO != 0) return INFO;

      LWORK = static_cast<int>(WORK[0]);
      LIWORK = IWORK[0];

      delete[] WORK;
      delete[] IWORK;

      WORK = new double[LWORK];
      IWORK = new int[LIWORK];
    }
    else if (backend == DSYGV) {
      LWORK = -1;
      WORK = new double[1];

      newscf_dsygv(&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB,
                   W, WORK, &LWORK, &INFO);
      if (INFO != 0) return INFO;

      LWORK = static_cast<int>(WORK[0]);
      delete[] WORK;
      WORK = new double[LWORK];
    }
    else if (backend == STDTRANS_SYEV) {
      // dpotrf and dtrsm do NOT require workspace arrays
      // Only dsyev requires workspace

      // Query dsyev workspace size
      int lwork_query = -1;
      double wkopt;
      newscf_dsyev(&JOBZ, &UPLO, &N, nullptr, &LDA, nullptr, &wkopt, &lwork_query, &INFO);
      if (INFO != 0) {
        std::cerr << "dsyev workspace query failed with INFO = " << INFO << std::endl;
        return INFO;
      }
      LWORK = static_cast<int>(wkopt);

      // Allocate one workspace array for dsyev
      WORK = new double[LWORK];
      // No IWORK needed here
      IWORK = nullptr;
      LIWORK = 0;
    }
    else {
      std::cerr << "Unknown backend!" << std::endl;
      return -1;
    }

    return INFO;
  }

  void set_A(const double* rA) {
    std::memcpy(A, rA, sizeof(double) * N * N);
  }

  void set_B(const double* rB) {
    std::memcpy(B, rB, sizeof(double) * N * N);
  }

  int run() {
    INFO = 0;

    if (backend == DSYGV) {
      newscf_dsygv(&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB,
                   W, WORK, &LWORK, &INFO);
      if (INFO != 0) {
        std::cerr << "DSYGV failed with INFO = " << INFO << std::endl;
        return INFO;
      }
    }
    else if (backend == DSYGVD) {
      newscf_dsygvd(&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB,
                    W, WORK, &LWORK, IWORK, &LIWORK, &INFO);
      if (INFO != 0) {
        std::cerr << "DSYGVD failed with INFO = " << INFO << std::endl;
        return INFO;
      }
    }
    else if (backend == STDTRANS_SYEV) {
      double* B_factor = new double[N * N];
      std::memcpy(B_factor, B, sizeof(double) * N * N);

      newscf_dpotrf(&UPLO, &N, B_factor, &LDB, &INFO);
      if (INFO != 0) {
        delete[] B_factor;
        return INFO;
      }

      char side = 'L';
      char diag = 'N';
      double alpha = 1.0;

      if (UPLO == 'U') {
        char trans = 'T';
        newscf_dtrsm(&side, &UPLO, &trans, &diag, &N, &N, &alpha, B_factor, &LDB, A, &LDA);
        trans = 'N';
        newscf_dtrsm(&side, &UPLO, &trans, &diag, &N, &N, &alpha, B_factor, &LDB, A, &LDA);
      } else {
        char trans = 'N';
        newscf_dtrsm(&side, &UPLO, &trans, &diag, &N, &N, &alpha, B_factor, &LDB, A, &LDA);
        trans = 'T';
        newscf_dtrsm(&side, &UPLO, &trans, &diag, &N, &N, &alpha, B_factor, &LDB, A, &LDA);
      }

      newscf_dsyev(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO);
      if (INFO != 0) {
        delete[] B_factor;
        return INFO;
      }

      if (JOBZ == 'V') {
        if (UPLO == 'U') {
          char trans = 'N';
          newscf_dtrsm(&side, &UPLO, &trans, &diag, &N, &N, &alpha, B_factor, &LDB, A, &LDA);
        } else {
          char trans = 'T';
          newscf_dtrsm(&side, &UPLO, &trans, &diag, &N, &N, &alpha, B_factor, &LDB, A, &LDA);
        }
      }

      delete[] B_factor;
    }

    else {
      std::cerr << "Unknown backend!" << std::endl;
      return -1;
    }

    return INFO;
  }

  void copy_eigenvectors(double* result) {
    std::memcpy(result, A, sizeof(double) * N * N);
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