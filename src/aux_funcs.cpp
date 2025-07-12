#include "newscf/ndtx.hpp"
#include "newscf/newscf.hpp"

namespace newscf {
	void build_rhf_density_matrix (IntegralEngineHandle* handle, ndtx::NDTX<double>& C, ndtx::NDTX<double>& D) {
		const int nbf = handle->nbf_tot;
		const int nocc = handle->nelec / 2;
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
		IntegralEngineHandle* handle,
		newscf::ndtx::NDTX<double>& H,
		newscf::ndtx::NDTX<double>& D,
		newscf::ndtx::NDTX<double>& F,
		newscf::ndtx::NDTX<double>& ERI) {
		const int nbf = handle->nbf_tot;
#pragma omp parallel for
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

	double energy_contraction (ndtx::NDTX<double>& D, ndtx::NDTX<double>& H, ndtx::NDTX<double>& F) {
		// TODO: rewrite using AVX+OpenMP
		double ener = 0;
		for (int i = 0; i < D.ndata; i++)
			ener += D.vectorGet(i) * (H.vectorGet(i) + F.vectorGet(i));
		return ener * 0.5;
	}

	// A * B * C -> D
	// A * B -> D
	// D * C -> D
	// A * B * C -> D, all square dim x dim
	void contract_matrix(ndtx::NDTX<double>& A, ndtx::NDTX<double>& B, ndtx::NDTX<double>& C, ndtx::NDTX<double>& D) {
		char TRANSA = 'N';
		char TRANSB = 'N';
		double ALPHA = 1.0;
		double BETA = 0.0;
		int dim = static_cast<int>(A.dims[0]);

		// Internal workspace
		ndtx::NDTX<double> W;
		W.resizeToMatrix(dim, dim);

		// W = A * B
		newscf_dgemm(&TRANSA, &TRANSB, &dim, &dim, &dim,
					 &ALPHA, A.data, &dim, B.data, &dim,
					 &BETA, W.data, &dim);

		// D = W * C
		newscf_dgemm(&TRANSA, &TRANSB, &dim, &dim, &dim,
					 &ALPHA, W.data, &dim, C.data, &dim,
					 &BETA, D.data, &dim);
	}


	// A B C - C B A
	void commutator_matrix (ndtx::NDTX<double> &A, ndtx::NDTX<double> &B, ndtx::NDTX<double> &C, ndtx::NDTX<double> &D) {
		// Internal Workspace
		ndtx::NDTX<double> W;
		W.resizeToMatrix(A.dims[0], A.dims[1]);

		contract_matrix(A, B, C, W);
		contract_matrix(C, B, A, D);

		// TODO: Parallelize and AVX
		for (int i = 0; i < A.dims[0] * A.dims[1]; i++) {
			D.data[i] = W.data[i] - D.data[i];
		}
	}
}