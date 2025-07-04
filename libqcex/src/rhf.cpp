#include <cstdlib>
#include <cstring>
#include <iomanip>

#include "qcex.hpp"
#include "qcex_utils.hpp"
#include "lapackinterface.hpp"
#include "logging_api.hpp"

extern "C" {
    #include "cint.h"
    #include "cint_funcs.h"
}

using newscf::libmatrix::Matrix;
using newscf::qcex::IntegralEngineHandle;

#ifndef LIBCINT_CACHE_SIZE
#define LIBCINT_CACHE_SIZE 2048
#endif

void _diagonalize_hcore (Matrix<double>* Hcore, Matrix<double>* overlap, Matrix<double>** coeffecients_matrix, Matrix<double>** eigenvalues) {
    
    int ITYPE = 1;
    char JOBZ = 'V';
    char UPLO = 'U';
    int N = Hcore -> sz_rows;
    double* A = new double [N * N];
    int LDA = Hcore -> sz_rows;
    double* B = new double [N * N];
    int LDB = overlap -> sz_cols;
    double* W = new double [N];
    double* WORK = new double [1];
    int LWORK = -1;
    int* IWORK = new int[1];
    int LIWORK = -1;
    int INFO = 0;

    memcpy (A, Hcore->data, sizeof(double) * N * N);
    memcpy (B, overlap->data, sizeof(double) * N * N);

    newscf_dsyevd (&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, IWORK, &LIWORK, &INFO);

    LWORK = static_cast<int>(WORK[0]);
    LIWORK = static_cast<int>(IWORK[0]);
    delete[] WORK, IWORK;
    WORK = new double[LWORK];
    IWORK = new int[LIWORK];

    newscf_dsyevd (&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, IWORK, &LIWORK, &INFO);
    
    if (INFO != 0)
        HANDLE_ERROR ("_diagonalize_hcore failed with INFO = "+std::to_string(INFO), 101);
    
    (*eigenvalues) = new Matrix<double>(W, N, 1);
    (*coeffecients_matrix) = new Matrix<double>(A, N, N);

    delete[] A, B, W, WORK, IWORK;

}

void _build_rhf_density_matrix (newscf::ndtx::NDTX<double>& C, newscf::ndtx::NDTX<double>& D, int nbf, int nocc) {
    for (int mu = 0; mu < nbf; ++mu) {
        for (int nu = 0; nu < nbf; ++nu) {
            double sum = 0.0;
            for (int i = 0; i < nocc; ++i) {
                sum += C.matrixGet(mu, i) * C.matrixGet(nu, i);
            }
            D.matrixSet(mu, nu, 2 * sum);
        }
    }
}

void _build_eri_tensor (IntegralEngineHandle* handle, newscf::ndtx::NDTX<double>& ERI) {
    const int nbf = handle->nbf_tot;
    const int nshells = handle->nbas;
    ERI.resizeToTensor4D(nbf);
    ERI.fill(0);

    CINTOpt* opt = nullptr;
    int2e_optimizer(&opt, handle->atm, handle->natm, handle->bas, handle->nbas, handle->env);
    double* cache = new double[LIBCINT_CACHE_SIZE];

    const int max_prim = handle->nbf_max * handle->nbf_max * handle->nbf_max * handle->nbf_max;
    double* buf = new double[max_prim];

    for (int i = 0; i < nshells; ++i) {
        int di = CINTcgto_spheric(i, handle->bas);
        int ioff = handle->shell_to_index[i];

        for (int j = 0; j <= i; ++j) {
            int dj = CINTcgto_spheric(j, handle->bas);
            int joff = handle->shell_to_index[j];

            for (int k = 0; k < nshells; ++k) {
                int dk = CINTcgto_spheric(k, handle->bas);
                int koff = handle->shell_to_index[k];

                for (int l = 0; l <= k; ++l) {
                    int dl = CINTcgto_spheric(l, handle->bas);
                    int loff = handle->shell_to_index[l];

                    // Enforce (ij) ≤ (kl) in lex order
                    if ((i * nshells + j) < (k * nshells + l)) continue;

                    int shls[4] = {i, j, k, l};
                    int dims[4] = {di, dj, dk, dl};

                    int ncomp = int2e_sph(buf, dims, shls,
                                          handle->atm, handle->natm,
                                          handle->bas, handle->nbas,
                                          handle->env, opt, cache);
                    if (ncomp != 1) HANDLE_ERROR("int2e_sph failed", 1);

                    for (int ii = 0; ii < di; ++ii) {
                        int mu = ioff + ii;
                        for (int jj = 0; jj < dj; ++jj) {
                            int nu = joff + jj;
                            for (int kk = 0; kk < dk; ++kk) {
                                int lam = koff + kk;
                                for (int ll = 0; ll < dl; ++ll) {
                                    int sig = loff + ll;

                                    int idx = (((ii * dj + jj) * dk + kk) * dl + ll);
                                    double val = buf[idx];

                                    ERI.tensor4DSet(mu, nu, lam, sig, val);
                                    ERI.tensor4DSet(nu, mu, lam, sig, val);
                                    ERI.tensor4DSet(mu, nu, sig, lam, val);
                                    ERI.tensor4DSet(nu, mu, sig, lam, val);
                                    ERI.tensor4DSet(lam, sig, mu, nu, val);
                                    ERI.tensor4DSet(lam, sig, nu, mu, val);
                                    ERI.tensor4DSet(sig, lam, mu, nu, val);
                                    ERI.tensor4DSet(sig, lam, nu, mu, val);

                                }
                            }
                        }
                    }
                }
            }
        }
    }

    CINTdel_optimizer(&opt);
    delete[] cache;
    delete[] buf;
}


void _compute_fock(
    IntegralEngineHandle* handle,
    newscf::ndtx::NDTX<double>& H,
    newscf::ndtx::NDTX<double>& D,
    newscf::ndtx::NDTX<double>& F
) {
    F = H;
    const int nshells = handle->nbas;
    const int nbf = handle->nbf_tot;

    CINTOpt* opt = nullptr;
    double* cache = new double[LIBCINT_CACHE_SIZE];

    cint2e_sph_optimizer(&opt, handle->atm, handle->natm, handle->bas, handle->nbas, handle->env);

    // Conservative allocation for max primitive size (e.g., s,p,d functions)
    const int max_nprim = handle->nbf_max * handle->nbf_max * handle->nbf_max * handle->nbf_max;
    double* buf_J = new double[max_nprim];
    double* buf_K = new double[max_nprim];

    for (int i = 0; i < nshells; ++i) {
        int di = CINTcgto_spheric(i, handle->bas);
        int ioff = handle->shell_to_index[i];

        for (int j = 0; j < nshells; ++j) {
            int dj = CINTcgto_spheric(j, handle->bas);
            int joff = handle->shell_to_index[j];

            for (int k = 0; k < nshells; ++k) {
                int dk = CINTcgto_spheric(k, handle->bas);
                int koff = handle->shell_to_index[k];

                for (int l = 0; l < nshells; ++l) {
                    int dl = CINTcgto_spheric(l, handle->bas);
                    int loff = handle->shell_to_index[l];

                    int nprim = di * dj * dk * dl;

                    int shls_J[4] = {i, j, k, l};
                    int dims_J[4] = {di, dj, dk, dl};

                    int ncomp_J = int2e_sph(buf_J, dims_J, shls_J,
                                            handle->atm, handle->natm,
                                            handle->bas, handle->nbas,
                                            handle->env, opt, cache);
                    if (ncomp_J != 1) HANDLE_ERROR("libcint 2e J failed", 102);

                    int shls_K[4] = {i, k, j, l};
                    int dims_K[4] = {di, dk, dj, dl};

                    int ncomp_K = int2e_sph(buf_K, dims_K, shls_K,
                                            handle->atm, handle->natm,
                                            handle->bas, handle->nbas,
                                            handle->env, opt, cache);
                    if (ncomp_K != 1) HANDLE_ERROR("libcint 2e K failed", 103);

                    for (int ii = 0; ii < di; ++ii) {
                        for (int jj = 0; jj < dj; ++jj) {
                            int mu = ioff + ii;
                            int nu = joff + jj;

                            double sum = 0.0;

                            for (int kk = 0; kk < dk; ++kk) {
                                for (int ll = 0; ll < dl; ++ll) {
                                    int la = koff + kk;
                                    int sig = loff + ll;

                                    double Dval = D.matrixGet(la, sig);
                                    if (Dval == 0.0) continue;

                                    int idx_J = (((ii * dj + jj) * dk + kk) * dl + ll);
                                    int idx_K = (((ii * dk + kk) * dj + jj) * dl + ll);

                                    double eri_J = buf_J[idx_J];
                                    double eri_K = buf_K[idx_K];

                                    sum += Dval * (eri_J - 0.5 * eri_K);
                                }
                            }

                            F.matrixSet(mu, nu, F.matrixGet(mu, nu) + sum);
                        }
                    }
                }
            }
        }
    }

    delete[] buf_J;
    delete[] buf_K;
    delete[] cache;
    CINTdel_optimizer(&opt);
}


namespace newscf::qcex {
    void  atom_rhf_energies   (IntegralEngineHandle* handle, double** energies) {
        //Matrix<double>* overlap = nullptr;
        /*Matrix<double>* nuclear = nullptr;
        Matrix<double>* kinetic = nullptr;
        calculate_overlap_matrix (handle, &overlap);
        calculate_kinetic_energy_matrix (handle, &kinetic);
        calculate_nuclear_attraction_matrix (handle, &nuclear);*/
        /*double* Tmat = nullptr;
        double* Smat = nullptr;
        calculate_hf_matrices(handle, &Tmat, &Smat);
        overlap = new Matrix<double>(Smat, handle->nbf_tot, handle->nbf_tot);
        Matrix<double> T = Matrix<double>(Tmat, handle->nbf_tot, handle->nbf_tot);
        Matrix<double>* C = nullptr;
        Matrix<double>* eval = nullptr;
        _diagonalize_hcore (&T, overlap, &C, &eval);
        Matrix<double> D (handle->nbf_tot, handle->nbf_tot);
        int nocc = 0;
        for (int i = 0; i < handle->natm; i++) nocc += handle->atm[6 * i]; // change this later fine for now
        std::cout << nocc << "\n";
        Matrix<double>* F = nullptr;
        // SCF loop
        std::cout << std::setprecision(16);
        for (int j = 0; j < 10; j++) {
            std::cout << "iter=" << (j+1);
            _build_rhf_density_matrix (C, &D, handle->nbf_tot, nocc);
            _compute_fock (handle, &T, &D, &F);
            _diagonalize_hcore (F, overlap, &C, &eval);
            for (int i = 0; i < handle->nbf_tot; i++) std::cout << ", " << eval->get(i, 0);
            std::cout << "\n";
        }*/

        NDTX<double> T; // Hcore
        NDTX<double> D; // Density Matrix
        NDTX<double> S; // Overlap Matrix
        NDTX<double> F; // Fock Tensor
        NDTX<double> C; // Coefficient Matrix

        _build_eri_tensor(handle, F);

        F.vectorPrint();

        /*int     ITYPE    = 1;
        char    JOBZ     = 'V';
        char    UPLO     = 'U';
        int     N        = handle->nbf_tot;
        int     LDA      = N;
        int     LDB      = N;
        double* W        = new double [N];
        double* WORK     = new double [1];
        int     LWORK    = -1;
        int     INFO     = 0;

        double* A = new double [N * N];
        double* B = new double [N * N];

        // Workspace Call
        newscf_dsygv (&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, &INFO);

        if (INFO != 0) HANDLE_ERROR("dsygv failed with INFO = "+std::to_string(INFO), 101);
        LWORK = static_cast<int>(WORK[0]);
        delete[] WORK;
        WORK = new double[LWORK];

        calculate_hf_matrices (handle, T, S);

        // Initial Hcore Guess
        std::copy (T.data, T.data + N * N, A);
        std::copy (S.data, S.data + N * N, B);
        newscf_dsygv (&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, &INFO);
        if (INFO != 0) HANDLE_ERROR("dsygv failed with INFO = "+std::to_string(INFO), 101);
        std::cout << "";
        for (int j = 0; j < N; ++j) std::cout << ", " << W[j];
        std::cout << std::endl;
        C.resizeToMatrix(N);
        std::copy (A, A + N * N, C.data);

        // Build Initial Fock Matrix
        D.resizeToMatrix(handle->nbf_tot);
        F.resizeToMatrix(handle->nbf_tot);
        _build_rhf_density_matrix (C, D, handle->nbf_tot, 2);

        for (int i = 0; i < 10; ++i) {
            _compute_fock (handle, T, D, F);
            std::copy (F.data, F.data + N * N, A);
            std::copy (S.data, S.data + N * N, B);
            newscf_dsygv (&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, &INFO);
            if (INFO != 0) HANDLE_ERROR("dsygv failed with INFO = "+std::to_string(INFO), 101);
            std::copy (A, A + N * N, C.data);
            std::cout << i;
            for (int j = 0; j < N; ++j) std::cout << ", " << W[j];
            std::cout << std::endl;
            _build_rhf_density_matrix (C, D, handle->nbf_tot, 2);
        }

        delete[] W;
        delete[] WORK;
        delete[] A;
        delete[] B;*/
    }
}