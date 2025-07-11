//
// Created by chlorinepentoxide on 10/7/25.
//

#ifndef DIIS_HPP
#define DIIS_HPP

#include "common.hpp"
#include "lapackinterface.hpp"
#include "ndtx.hpp"

#include <cfloat>

namespace newscf {
	enum DIIS_TYPE {
		CDIIS, EDIIS, ADIIS
	};

	class DIIS {
	public:
		double* D_data;  // Density Matrix Data
		double* F_data;  // Fock Matrix Data
		double* S_data;  // Overlap Matrix Data
		double* E_data;  // Energy Data
		int     dim;     // Matrix Dimension
		int     nhist;   // Number of previous matrices to store
		int     cptr;    // (Current) Matrix Loc Pointer; resets when to when nhist is hit
		int     gcptr;   // Global Counter; does not reset
		DIIS_TYPE type;
		GeneralLinsolve_Backend backend;

		// LAPACK/BLAS and Internal Workspace Parameters
		double* coeff_vector; // output of any of run_xDIIS
		double* work_i1; // BLAS matrix i (left commutator); intermediate in B_ij = <work_i|work_j>
		double* work_i2; // BLAS matrix i (right commutator); intermediate in B_ij = <work_i|work_j>
		double* work_j1; // BLAS matrix j (left commutator); intermediate in B_ij = <work_i|work_j>
		double* work_j2; // BLAS matrix j (right commutator); intermediate in B_ij = <work_i|work_j>
		ndtx::NDTX<double> B; // Residual Matrix

		// BLAS dgemm call parameters
		char TRANSA = 'N';
		char TRANSB = 'N';
		double ALPHA = 1.0;
		double BETA  = 0.0;

		DIIS() {
			D_data = nullptr;
			F_data = nullptr;
			S_data = nullptr;
			E_data = nullptr;
			coeff_vector = nullptr;
			work_i1 = nullptr;
			work_i2 = nullptr;
			work_j1 = nullptr;
			work_j2 = nullptr;
			dim = 1;
			nhist = 1;
			cptr = 0;
			gcptr = 0;
			type = ADIIS;
		}

		~DIIS() {
			delete[] D_data;
			delete[] F_data;
			delete[] E_data;
			delete[] work_i1;
			delete[] work_i2;
			delete[] work_j1;
			delete[] work_j2;
		}

		inline void set_backend (GeneralLinsolve_Backend b) {
			backend = b;
		}

		inline void set_diis_type(const DIIS_TYPE t) {
			type = t;
		}

		inline void set_diis_dim(const int d) {
			dim = d;
		}

		inline void set_diis_nhist(const int n) {
			nhist = n;
		}

		inline void init() {
			coeff_vector = new double[nhist + 1];
			memset(coeff_vector, 0, nhist * sizeof(double));
			work_i1 = new double[dim * dim];
			work_i2 = new double[dim * dim];
			work_j1 = new double[dim * dim];
			work_j2 = new double[dim * dim];
			D_data = new double[nhist * dim * dim];
			memset(D_data, 0, nhist * dim * dim * sizeof(double));
			F_data = new double[nhist * dim * dim];
			memset(F_data, 0, nhist * dim * dim * sizeof(double));
			E_data = new double[nhist];
			memset(E_data, 0, nhist * sizeof(double));
		}

		inline void set_S (const ndtx::NDTX<double>& rS) {
			S_data = rS.data;
		}

		inline void store (const ndtx::NDTX<double>& rD, const ndtx::NDTX<double>& rF, const double E) {
			memcpy(D_data + cptr * dim * dim, rD.data, dim * dim * sizeof(double));
			memcpy(F_data + cptr * dim * dim, rF.data, dim * dim * sizeof(double));

			E_data[cptr] = E;

			cptr = (cptr + 1) % nhist;
			gcptr++;
		}

		inline void run_CDIIS() {
			int lim = (gcptr >= nhist) ? nhist : gcptr;
			B.resizeToMatrix(lim);

			// To compute the B matrix
			for (int i = 0; i < lim; i++) {
				for (int j = 0; j <= i; j++) {

					// ensure clean workspace
					memset(work_i1, 0, dim * dim * sizeof(double));
					memset(work_i2, 0, dim * dim * sizeof(double));
					memset(work_j1, 0, dim * dim * sizeof(double));
					memset(work_j2, 0, dim * dim * sizeof(double));

					double* D_mat_i = D_data + i * dim * dim;
					double* F_mat_i = F_data + i * dim * dim;
					double* D_mat_j = D_data + j * dim * dim;
					double* F_mat_j = F_data + j * dim * dim;

					// F_i * D_i -> work_i1
					newscf_dgemm(&TRANSA, &TRANSB, &dim, &dim, &dim, &ALPHA, F_mat_i, &dim, D_mat_i, &dim, &BETA, work_i1, &dim);
					// D_i * F_i -> work_i2
					newscf_dgemm(&TRANSA, &TRANSB, &dim, &dim, &dim, &ALPHA, D_mat_i, &dim, F_mat_i, &dim, &BETA, work_i2, &dim);
					// work_i1/FD * S -> work_i1
					newscf_dgemm(&TRANSA, &TRANSB, &dim, &dim, &dim, &ALPHA, work_i1, &dim, S_data, &dim, &BETA, work_i1, &dim);
					// S * work_i2/DF -> work_i2
					newscf_dgemm(&TRANSA, &TRANSB, &dim, &dim, &dim, &ALPHA, S_data, &dim, work_i2, &dim, &BETA, work_i2, &dim);

					// F_j * D_j -> work_j1
					newscf_dgemm(&TRANSA, &TRANSB, &dim, &dim, &dim, &ALPHA, F_mat_j, &dim, D_mat_j, &dim, &BETA, work_j1, &dim);
					// D_j * F_j -> work_j2
					newscf_dgemm(&TRANSA, &TRANSB, &dim, &dim, &dim, &ALPHA, D_mat_j, &dim, F_mat_j, &dim, &BETA, work_j2, &dim);
					// work_j1/FD * S -> work_j1
					newscf_dgemm(&TRANSA, &TRANSB, &dim, &dim, &dim, &ALPHA, work_j1, &dim, S_data, &dim, &BETA, work_j1, &dim);
					// S * work_j2/DF -> work_j2
					newscf_dgemm(&TRANSA, &TRANSB, &dim, &dim, &dim, &ALPHA, S_data, &dim, work_j2, &dim, &BETA, work_j2, &dim);

					// Contract work_i1 - work_i2 (FDS - SDF)_i -> work_i1
					for (int k = 0; k < dim * dim; k++) work_i1[k] -= work_i2[k];
					// Contract work_j1 - work_j2 (FDS - SDF)_j -> work_j1
					for (int k = 0; k < dim * dim; k++) work_j1[k] -= work_j2[k];

					// Dot
					double dot_prod = 0.0;
					for (int k = 0; k < dim * dim; k++) dot_prod += work_i1[k] * work_j1[k];

					// B_ij = B_ji
					B.matrixSet(i, j, dot_prod);
					B.matrixSet(j, i, dot_prod);
				}
			}

			// Holy Grail of Linear Algebra; Solve X for AX=B
			ndtx::NDTX<double> linsolve_A; linsolve_A.resizeToMatrix(lim + 1);
			ndtx::NDTX<double> linsolve_B; linsolve_B.resizeToVector(lim + 1);

			// Embed residual matrix B inside linsolve_A
			for (int i = 0; i < lim; i++) for (int j = 0; j < lim; j++) linsolve_A.matrixSet(i, j, B.matrixGet(i, j));

			// Decorate matrix linsolve_A & linsolve_B
			for (int i = 0; i < lim; i++) {
				linsolve_A.matrixSet(i, lim, -1);
				linsolve_A.matrixSet(lim, i, -1);
				linsolve_B.vectorSet(i, 0);
			}

			linsolve_A.matrixSet(lim, lim, 0);
			linsolve_B.vectorSet(lim, 1);

			// Solve
			GeneralLinsolve linsolver;
			linsolver.set_N (lim+1);
			linsolver.set_NRHS (1);
			linsolver.set_backend (DSYSV);
			linsolver.init();

			linsolver.set_A(linsolve_A.data);
			linsolver.set_B(linsolve_B.data);
			linsolver.run();
			linsolver.copy_result(coeff_vector);
		}

		inline void run_EDIIS() {
			const double delta = 1e-8;
			double minE = DBL_MAX;
			const int lim = (gcptr >= nhist) ? nhist : gcptr;

			// Find minimum energy
			for (int i = 0; i < lim; i++) {
				if (E_data[i] < minE) minE = E_data[i];
			}

			// Compute weights: inverse of shifted energies
			double sum = 0.0;
			for (int i = 0; i < lim; i++) {
				coeff_vector[i] = 1.0 / (E_data[i] - minE + delta);
				sum += coeff_vector[i];
			}

			// Normalize coefficients so sum to 1 (convex combination)
			for (int i = 0; i < lim; i++) {
				coeff_vector[i] /= sum;
			}
		}


		inline void build_next (ndtx::NDTX<double>& nF) {
			const int lim = (gcptr >= nhist) ? nhist : gcptr;
			for (int i = 0; i < dim * dim; i++) {
				double val = 0;
				for (int k = 0; k < lim; k++) {
					double coeff = coeff_vector[k];
					int idx = k * dim * dim + i;
					val += coeff * F_data[idx];
				}
				nF.vectorSet(i, val);
			}
		}
	};
}

#endif //DIIS_HPP
