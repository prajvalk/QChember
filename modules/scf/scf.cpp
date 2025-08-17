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
#include <cfloat>

#include "newscf/newscf.hpp"
#include "newscf/cis.hpp"

namespace newscf {

    void build_rhf_density_matrix (const IntegralEngineHandle<LIBCINT>& handle, LACIS<DENSE_LAPACK, double>& C, LACIS<DENSE_LAPACK, double>& D) {
        const int nbf = handle.nbf_tot;
        const int nocc = handle.nelec / 2;
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

    void build_rhf_fock_matrix (
        const IntegralEngineHandle<LIBCINT>& handle,
        LACIS<DENSE_LAPACK, double>& H,
        LACIS<DENSE_LAPACK, double>& D,
        LACIS<DENSE_LAPACK, double>& F,
        LACIS<DENSE_LAPACK, double>& ERI) {
        const int nbf = handle.nbf_tot;
        for (int lam = 0; lam < nbf; lam++) {
            for (int sig = 0; sig <= lam; ++sig) {
                double val = 0;
                for (int mu = 0; mu < nbf; ++mu) {
                    for (int nu = 0; nu < nbf; ++nu) {
                        double J = ERI.tensor4DGet(mu, nu, lam, sig);
                        double K = ERI.tensor4DGet(mu, lam, nu, sig);
                        val += D.matrixGet(mu, nu) * (J - 0.5 * K);
                    }
                }
                F.matrixSet(lam, sig, H.matrixGet(lam, sig) + val);
                F.matrixSet(sig, lam, H.matrixGet(lam, sig) + val);
            }
        }
    }

    double energy_contraction (LACIS<DENSE_LAPACK, double>& D, LACIS<DENSE_LAPACK, double>& H, LACIS<DENSE_LAPACK, double>& F) {
        double ener = 0;
        for (int i = 0; i < D.ndata; i++)
            ener += D.vectorGet(i) * (H.vectorGet(i) + F.vectorGet(i));
        return ener * 0.5;
    }

    template<>
    double run_rhf<LIBCINT, DENSE_LAPACK>(const IntegralEngineHandle<LIBCINT>& handle) {
        LACIS<DENSE_LAPACK, double> T; // H_core Matrix
        LACIS<DENSE_LAPACK, double> S; // Overlap Matrix
        LACIS<DENSE_LAPACK, double> D; // Density Matrix
        LACIS<DENSE_LAPACK, double> ERI; // ERI Tensor
        LACIS<DENSE_LAPACK, double> F; // Fock Matrix
        LACIS<DENSE_LAPACK, double> C; // Coefficient Matrix

        T.resizeToMatrix(handle.nbf_tot);
        S.resizeToMatrix(handle.nbf_tot);
        D.resizeToMatrix(handle.nbf_tot);
        ERI.resizeToTensor4D(handle.nbf_tot);
        F.resizeToMatrix(handle.nbf_tot);
        C.resizeToMatrix(handle.nbf_tot);

        GeneralizedEigenvalueSolver<DSYGV> eigensolver;
        eigensolver.I1 = 1;
        eigensolver.CH1 = 'V';
        eigensolver.CH2 = 'U';
        eigensolver.set_n(handle.nbf_tot);
        eigensolver.init();

        // Initial Guess
        cint1e_hf(handle, T, S);
        cint2e_hf_eri(handle, ERI);
        eigensolver.set_a(T);
        eigensolver.set_b(S);
        eigensolver.solve();
        eigensolver.copy_eigenvectors(C);

        build_rhf_density_matrix(handle, C, D);

        double prev_ener = -DBL_MAX;
        bool has_converged = false;

        while (!has_converged) {
            build_rhf_fock_matrix(handle, T, D, F, ERI);
            eigensolver.set_a(F);
            eigensolver.set_b(S);
            eigensolver.solve();
            eigensolver.copy_eigenvectors(C);
            build_rhf_density_matrix(handle, C, D);
            double ener = energy_contraction(D, T, F);
            std::cout << "ener: " << ener << std::endl;
            if (fabs(prev_ener - ener) < 1e-10) has_converged = true;
            prev_ener = ener;
        }

        return prev_ener;
    }

}