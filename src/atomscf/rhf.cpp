#include "newscf/newscf.hpp"
#include "newscf/ndtx.hpp"
#include "newscf/common.hpp"
#include "newscf/lapackinterface.hpp"
#include "newscf/diis.hpp"

#include <cmath>
#include <cfloat>

namespace newscf {

	void atomscf_rhf (const SolverOptions opts, HFResult* result, IntegralEngineHandle* handle) {
		using ndtx::NDTX;

		NDTX<double> T;  // Hcore
		NDTX<double> D;  // Density Matrix
		NDTX<double> S;  // Overlap Matrix
		NDTX<double> ERI;// ERI Tensor
		NDTX<double> F;  // Fock Matrix
		NDTX<double> C;  // Coefficient Matrix

		// DIIS
		DIIS diis;
		diis.set_diis_nhist(5);
		diis.set_diis_type(CDIIS);
		diis.set_diis_dim(handle->nbf_tot);
		diis.init();

		GeneralGEV eigensolver; // General eigensolver
		eigensolver.set_backend(opts.HF_GEV_BACKEND);

		calculate_eri_tensor (handle, ERI);
		calculate_hf_matrices (handle, T, S);

		diis.set_S(S);

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

		bool converged = false;
		double oldE = 0, rmsD = DBL_MAX;
		for (int i = 0; i < opts.HF_MAX_ITER; i++) {
			// Build Fock Matrix
			build_rhf_fock_matrix(handle, T, D, F, ERI);

			diis.store (D, F);

			if (i >= 3 && rmsD >= 1e-2) {
				diis.run_CDIIS();
				diis.build_next(F);
			}

			// Eigendecomposition
			eigensolver.set_A(F.data);
			eigensolver.set_B(S.data);
			eigensolver.run();
			eigensolver.copy_eigenvectors(C.data);

			// Rebuild Density Matrix
			const NDTX<double> D_old = D;
			build_rhf_density_matrix(handle, C, D);

			// Calculate Energy
			double ener = energy_contraction(D, T, F);

			// Check Convergence
			rmsD = D_old.rms(D);
			double deltaE = fabs(ener - oldE);
			oldE = ener;

			if (opts.HF_STRONG_CONV && rmsD < opts.HF_CONV_DTOL && deltaE < opts.HF_CONV_DTOL) {
				converged = true;
				break;
			} else if (!opts.HF_STRONG_CONV && (rmsD < opts.HF_CONV_DTOL || deltaE < opts.HF_CONV_DTOL)) {
				converged = true;
				break;
			}

			LOG (INFO, "iter="+std::to_string(i)+" ener="+std::to_string(ener)+" delE="+std::to_string(deltaE)+" rmsD="+std::to_string(rmsD));
			if (opts.VERBOSE) std::cout << "iter="+std::to_string(i)+" ener="+std::to_string(ener)+" delE="+std::to_string(deltaE)+" rmsD="+std::to_string(rmsD) << std::endl;
		}

		if (!converged) {
			LOG (WARN, "RHF SCF procedure has not converged within HF_MAX_ITER");
			LOG (WARN, "Density Matrix RMS Error = "+std::to_string(rmsD));
		}

		result->HF_ENER = oldE;
	}
}