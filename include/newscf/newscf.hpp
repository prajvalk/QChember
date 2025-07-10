/**
 * NewSCF Library API Header
 * Mozilla Public License 2.0
 * Copyright (C) 2025, Prajval K
 */

#ifndef NEWSCF_HPP
#define NEWSCF_HPP

#include "ndtx.hpp"


#include <string>

namespace newscf {

	int load_molecule (const std::string& fn, double** space);

	void destroy_molecule (double* molecule);

	// Loads a complete basis from file to memory; Only Psi4 format is supported
	int  load_basis       (const std::string& filename, double** basis);
	// Loads a partial basis; only basis functions relating to the molecule information
	int  load_basis       (std::string filename, double** basis, double* molecule);
	// Destroys all memory associated with the basis
	void destroy_basis    (double* basis);

	struct IntegralEngineHandle {
		int*    atm;       // Atomic Description
		int*    bas;       // Basis Description
		double* env;       // Environment Description
		int*    shell_to_index; // Internal Shell Indexing Map
		double* cache;     // libcint cache
		int     natm;      // Number of Atoms
		int     nbas;      // Number of Basis Shells
		int     nenv;
		int     nelec;     // Total number of electrons in system
		int     nbf_tot;   // Number of Basis Functions (Total)
		int     nbf_max;   // Number of Basis Functions (Max in any shell)
	};

	void  initialize_handle (double* molecule, int molecule_size, double* basis, int basis_size, IntegralEngineHandle* handle);
	void  destroy_handle   (IntegralEngineHandle* handle);

	void calculate_hf_matrices                 (IntegralEngineHandle* handle, ndtx::NDTX<double>& T, ndtx::NDTX<double>& S);
	void calculate_eri_tensor                  (IntegralEngineHandle* handle, ndtx::NDTX<double>& ERI);

	void build_rhf_density_matrix              (IntegralEngineHandle* handle, ndtx::NDTX<double>& C, ndtx::NDTX<double>& D);
	void build_rhf_fock_matrix				   (IntegralEngineHandle* handle, ndtx::NDTX<double>& H, ndtx::NDTX<double>& D, ndtx::NDTX<double>& F, ndtx::NDTX<double>& ERI);

	// ATOM SCF

	struct SolverOptions {
		double HF_CONV_ETOL = 1e-10;
		double HF_CONV_DTOL = 1e-10;
		double HF_MAX_ITER  = 100;
		bool   HF_STRONG_CONV = true;
	};

	struct HFResult {
		double HF_ENER = 0;
	};

	void atomscf_rhf (SolverOptions opts, HFResult* result, IntegralEngineHandle* handle);

}

#endif //NEWSCF_HPP