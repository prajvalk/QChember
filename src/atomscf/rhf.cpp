#include "newscf/newscf.hpp"
#include "newscf/ndtx.hpp"
#include "newscf/common.hpp"
#include "newscf/lapackinterface.hpp"

#include <cmath>

namespace newscf {

	void atomscf_rhf (const SolverOptions opts, IntegralEngineHandle* handle) {
		using ndtx::NDTX;

		NDTX<double> T;  // Hcore
		NDTX<double> D;  // Density Matrix
		NDTX<double> S;  // Overlap Matrix
		NDTX<double> ERI;// ERI Tensor
		NDTX<double> F;  // Fock Matrix
		NDTX<double> C;  // Coefficient Matrix

		calculate_eri_tensor (handle, ERI);
		calculate_hf_matrices (handle, T, S);

		if (handle->nelec % 2 == 1) {
			LOG (WARN, "You are trying to run a RHF calculation with an odd electron count (nelec="+std::to_string(handle->nelec)+"). Consider using ROHF or UHF instead.");
		}

		int nocc = handle->nelec / 2;

		// Workspace Preparation

		int     ITYPE    = 1;
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

		// Initial Hcore Guess
        std::copy (T.data, T.data + N * N, A);
        std::copy (S.data, S.data + N * N, B);

        newscf_dsygv (&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, &INFO);
		S.matrixPrint();
        if (INFO != 0) HANDLE_ERROR("dsygv failed with INFO = "+std::to_string(INFO), 101);

        std::cout << "";
        for (int j = 0; j < N; ++j) std::cout << ", " << W[j];
        std::cout << std::endl;
        C.resizeToMatrix(N);
        std::copy (A, A + N * N, C.data);

        // Build Initial Fock Matrix
        D.resizeToMatrix(N);
        F.resizeToMatrix(N);
        build_rhf_density_matrix (handle, C, D);

        double energy = 0.0;
        double last_energy = 0.0;
        double rms_D = 1.0;
        int max_iter = 100;
        int iter = 0;
        double conv_thresh = 1e-12;

        while (iter < max_iter && rms_D > conv_thresh) {
            // Build Fock matrix
            build_rhf_fock_matrix(handle, T, D, F, ERI);

            // Solve F C = S C E
            std::copy(F.data, F.data + N*N, A);
            std::copy(S.data, S.data + N*N, B);
            newscf_dsygv(&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, &INFO);
            if (INFO != 0) HANDLE_ERROR("dsygv failed", 1);

            std::copy(A, A + N*N, C.data);

            // Save old D
            NDTX<double> D_old = D;

            // Build new D
            build_rhf_density_matrix(handle, C, D);

            // Compute electronic energy
            double E_elec = 0.0;
            for (int mu = 0; mu < N; ++mu) {
                for (int nu = 0; nu < N; ++nu) {
                    E_elec += D.matrixGet(mu, nu) * (T.matrixGet(mu, nu) + F.matrixGet(mu, nu));
                }
            }
            energy = E_elec * 0.5;

            // Compute RMS of D
            rms_D = 0.0;
            for (int i = 0; i < N * N; ++i) {
                double delta = D.data[i] - D_old.data[i];
                rms_D += delta * delta;
            }
            rms_D = std::sqrt(rms_D);

            std::cout << "iter: " << iter
                      << ", energy: " << energy
                      << ", deltaE: " << std::abs(energy - last_energy)
                      << ", rmsD: " << rms_D << std::endl;

            last_energy = energy;
            ++iter;
        }

		delete[] WORK;
		delete[] A;
		delete[] B;
		delete[] W;
	}
}