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
}