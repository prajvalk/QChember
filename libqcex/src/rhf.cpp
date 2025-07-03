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

void _build_rhf_density_matrix (Matrix<double>* coeff_matrix, Matrix<double>* density_matrix, int nbf, int nocc) {
    for (int mu = 0; mu < nbf; ++mu) {
        for (int nu = 0; nu < nbf; ++nu) {
            double sum = 0.0;
            for (int i = 0; i < nocc; ++i) {
                sum += coeff_matrix->get(mu, i) * coeff_matrix->get(nu, i);
            }
            density_matrix->set(mu, nu, 2 * sum);
        }
    }
}

void _compute_fock (IntegralEngineHandle* handle, Matrix<double>* Hcore, Matrix<double>* D, Matrix<double>** F) {
    (*F) = new Matrix<double>(handle->nbf_tot, handle->nbf_tot);
    const int nshells = handle->nbas;

    CINTOpt* opt = nullptr;
    double* cache = new double[LIBCINT_CACHE_SIZE];

    cint2e_sph_optimizer(&opt, handle->atm, handle->natm, handle->bas, handle->nbas, handle->env);

    for (int i = 0; i < nshells; ++i) {
        int di = CINTcgto_spheric(i, handle->bas);
        int ioff = handle->shell_to_index[i];
        for (int j = 0; j <= i; ++j) { // exploit symmetry
            int dj = CINTcgto_spheric(j, handle->bas);
            int joff = handle->shell_to_index[j];
            for (int k = 0; k < nshells; ++k) {
                int dk = CINTcgto_spheric(k, handle->bas);
                int koff = handle->shell_to_index[k];
                for (int l = 0; l < nshells; ++l) {
                    int dl = CINTcgto_spheric(l, handle->bas);
                    int loff = handle->shell_to_index[l];

                    int nprim = di * dj * dk * dl;
                    double* buf = new double[nprim];

                    int shls[4] = {i, j, k, l};
                    int dims[4] = {di, dj, dk, dl};

                    int ncomp = int2e_sph(buf, dims, shls,
                                        handle->atm, handle->natm,
                                        handle->bas, handle->nbas,
                                        handle->env, opt, cache);
                    
                    if (ncomp != 1) HANDLE_ERROR ("libcint 2e computation failed with ncomp = "+std::to_string(ncomp), 102);

                    for (int ii = 0; ii < di; ++ii) {
                        for (int jj = 0; jj < dj; ++jj) {
                            int mu = ioff + ii;
                            int nu = joff + jj;
                            double sum = 0.0;

                            for (int kk = 0; kk < dk; ++kk) {
                                for (int ll = 0; ll < dl; ++ll) {
                                    int la = koff + kk;
                                    int sig = loff + ll;

                                    int idx = (((ii * dj + jj) * dk + kk) * dl + ll);
                                    double eri_ij_kl = buf[idx];

                                    // For exchange, reorder (μλ|νσ) = (μν|λσ) with μ ↔ λ
                                    int idx_K = (((ii * dk + kk) * dj + jj) * dl + ll);
                                    double eri_il_kj = buf[idx_K];

                                    double Dval = D->get(la, sig);
                                    sum += Dval * (eri_ij_kl - 0.5 * eri_il_kj);
                                }
                            }

                            double Fval = Hcore->get(mu, nu) + sum;
                            (*F)->set(mu, nu, Fval);
                            if (mu != nu)
                                (*F)->set(nu, mu, Fval);  // Symmetrize
                        }
                    }
                    delete[] buf;
                }
            }
        }
    }

    delete[] cache;
    CINTdel_optimizer(&opt);

}

namespace newscf::qcex {
    void  atom_rhf_energies   (IntegralEngineHandle* handle, double** energies) {
        Matrix<double>* overlap = nullptr;
        Matrix<double>* nuclear = nullptr;
        Matrix<double>* kinetic = nullptr;
        calculate_overlap_matrix (handle, &overlap);
        calculate_kinetic_energy_matrix (handle, &kinetic);
        calculate_nuclear_attraction_matrix (handle, &nuclear);
        Matrix<double> T = (*nuclear) + (*kinetic);
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
        }
    }
}