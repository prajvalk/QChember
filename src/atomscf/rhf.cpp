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

		DSYGV eigensolver; // DYSGV eigensolver

		calculate_eri_tensor (handle, ERI);
		calculate_hf_matrices (handle, T, S);

		if (handle->nelec % 2 == 1) {
			LOG (WARN, "You are trying to run a RHF calculation with an odd electron count (nelec="+std::to_string(handle->nelec)+"). Consider using ROHF or UHF instead.");
		}

		// Workspace Preparation
		eigensolver.compute_eigenvectors(true);
		eigensolver.set_N(handle->nbf_tot);
		eigensolver.init();
		C.resizeToMatrix(handle->nbf_tot);
		D.resizeToMatrix(handle->nbf_tot);
		F.resizeToMatrix(handle->nbf_tot);

		// Initial Hcore Guess
		eigensolver.set_A(T.data);
		eigensolver.set_B(S.data);
		eigensolver.run();
		eigensolver.copy_eigenvectors(C.data);

		// Build Initial Density Matrix
		build_rhf_density_matrix(handle, C, D);

		for (int i = 0; i < 100; i++) {
			// Build Fock Matrix
			build_rhf_fock_matrix(handle, T, D, F, ERI);

			// Eigendecomposition
			eigensolver.set_A(F.data);
			eigensolver.set_B(S.data);
			eigensolver.run();
			eigensolver.copy_eigenvectors(C.data);

			// Rebuild Density Matrix
			const NDTX<double> D_old = D;
			build_rhf_density_matrix(handle, C, D);

			// Calculate Energy
			double ener = 0;
			for (int j = 0; j < handle->nbf_tot * handle->nbf_tot; j++) {
				ener += D.vectorGet(j) * (T.vectorGet(j) + F.vectorGet(j));
			}
			ener *= 0.5;

			std::cout << ener << std::endl;
		}
	}
}