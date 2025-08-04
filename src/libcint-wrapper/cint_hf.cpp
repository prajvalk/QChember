/*
 * NewSCF
 * Copyright (C) 2025, Prajval K
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cmath>

extern "C" {
#include "cint.h"
#include "cint_funcs.h"
}

#include "newscf/cis.hpp"
#include "newscf/newscf.hpp"

namespace newscf {

    template<>
    void cint1e_hf(const IntegralEngineHandle<LIBCINT>& handle, LACIS<DENSE_LAPACK, double>& T, LACIS<DENSE_LAPACK, double>& S) {
        using namespace newscf;
        double* buffer_K = new double[handle.nbf_max * handle.nbf_max];
        double* buffer_N = new double[handle.nbf_max * handle.nbf_max];
        double* buffer_S = new double[handle.nbf_max * handle.nbf_max];
        int*    shells = new int   [2];
        int*    dims   = new int   [3] {1, 1, 1};
        T.resizeToMatrix(handle.nbf_tot);
        T.fill(0);
        S.resizeToMatrix(handle.nbf_tot);
        T.fill(0);
        CINTOpt* opt_K = nullptr;
        CINTOpt* opt_N = nullptr;
        CINTOpt* opt_S = nullptr;
        int1e_kin_optimizer (&opt_K, handle.atm, handle.natm, handle.bas, handle.nbas, handle.env);
        int1e_nuc_optimizer (&opt_N, handle.atm, handle.natm, handle.bas, handle.nbas, handle.env);
        int1e_ovlp_optimizer(&opt_S, handle.atm, handle.natm, handle.bas, handle.nbas, handle.env);
        for (int i = 0; i < handle.nbas; i++) {
            for (int j = 0; j <= i; j++) {
                shells[0] = i;
                shells[1] = j;
                dims[0]   = CINTcgto_spheric(i, handle.bas);
                dims[1]   = CINTcgto_spheric(j, handle.bas);
                int ioff  = handle.shell_to_index[i];
                int joff  = handle.shell_to_index[j];
                memset(buffer_K, 0, sizeof(double) * dims[0] * dims[1]);
                memset(buffer_N, 0, sizeof(double) * dims[0] * dims[1]);
                memset(buffer_S, 0, sizeof(double) * dims[0] * dims[1]);
                int ncomp_K = int1e_kin_sph  (buffer_K, dims, shells, handle.atm, handle.natm, handle.bas, handle.nbas, handle.env, opt_K, handle.cache);
                int ncomp_N = int1e_nuc_sph  (buffer_N, dims, shells, handle.atm, handle.natm, handle.bas, handle.nbas, handle.env, opt_N, handle.cache);
                int ncomp_S = int1e_ovlp_sph (buffer_S, dims, shells, handle.atm, handle.natm, handle.bas, handle.nbas, handle.env, opt_S, handle.cache);
                if (ncomp_K != 1 || ncomp_N != 1 || ncomp_S != 1) HANDLE_ERROR ("libcint 1e calculation failed.", 101);
                for (int ii = 0; ii < dims[0]; ii++) {
                    for (int jj = 0; jj < dims[1]; jj++) {
                        double Kval = buffer_K[ii * dims[1] + jj];
                        double Nval = buffer_N[ii * dims[1] + jj];
                        double Sval = buffer_S[ii * dims[1] + jj];
                        T.matrixSet(ioff + ii, joff + jj, Kval + Nval);
                        T.matrixSet(joff + jj, ioff + ii, Kval + Nval);
                        S.matrixSet(ioff + ii, joff + jj, Sval);
                        S.matrixSet(joff + jj, ioff + ii, Sval);
                    }
                }
            }
        }
        CINTdel_optimizer (&opt_K);
        CINTdel_optimizer (&opt_N);
        CINTdel_optimizer (&opt_S);
        delete[] buffer_K;
        delete[] buffer_N;
        delete[] buffer_S;
        delete[] dims;
        delete[] shells;
    }

    template<>
    void cint2e_hf_eri<LIBCINT, DENSE_LAPACK, double> (const IntegralEngineHandle<LIBCINT>& handle, LACIS<DENSE_LAPACK, double>& ERI) {
	    // TODO: Parallelize via OpenMP

        const int nbf = handle.nbf_tot;
        const int nshells = handle.nbas;
        ERI.resizeToTensor4D(nbf);
        ERI.fill(0);

        CINTOpt* opt = nullptr;
        int2e_optimizer(&opt, handle.atm, handle.natm, handle.bas, handle.nbas, handle.env);
        double* cache = handle.cache;

        const int max_prim = handle.nbf_max * handle.nbf_max * handle.nbf_max * handle.nbf_max;
        double* buf = new double[max_prim];

        for (int i = 0; i < nshells; ++i) {
            int di = CINTcgto_spheric(i, handle.bas);
            int ioff = handle.shell_to_index[i];

            for (int j = 0; j <= i; ++j) {
                int dj = CINTcgto_spheric(j, handle.bas);
                int joff = handle.shell_to_index[j];

                for (int k = 0; k < nshells; ++k) {
                    int dk = CINTcgto_spheric(k, handle.bas);
                    int koff = handle.shell_to_index[k];

                    for (int l = 0; l <= k; ++l) {
                        int dl = CINTcgto_spheric(l, handle.bas);
                        int loff = handle.shell_to_index[l];

                        // Enforce (ij) ≤ (kl) in lex order
                        if ((i * nshells + j) < (k * nshells + l)) continue;

                        int shls[4] = {i, j, k, l};
                        int dims[4] = {di, dj, dk, dl};

                        int ncomp = int2e_sph(buf, dims, shls,
                                              handle.atm, handle.natm,
                                              handle.bas, handle.nbas,
                                              handle.env, opt, cache);
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

        CINTdel_2e_optimizer(&opt);
        delete[] buf;
    }

}