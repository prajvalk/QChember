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

#ifndef NEWSCF_HPP
#define NEWSCF_HPP

#include <string>
#include "cis.hpp"

namespace newscf {
	using namespace newscf::cis;

    int load_molecule     (const std::string& fn, double** space);

	void destroy_molecule (double* molecule);

	// Loads a complete basis from file to memory; Only Psi4 format is supported
	int  load_basis       (const std::string& filename, double** basis);

	// Loads a partial basis; only basis functions relating to the molecule information
	int  load_basis       (const std::string& filename, double** basis, const double* geom, int num_atoms);

	// Destroys all memory associated with the basis
	void destroy_basis    (double* basis);

	enum IntegralEngine { LIBCINT };

	template<IntegralEngine>
	class IntegralEngineHandle {
	public:
		void*   W1;
		void*   W2;
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

		IntegralEngineHandle() {
			W1  = nullptr;
			W2  = nullptr;
			atm = nullptr;
			bas = nullptr;
			env = nullptr;
			shell_to_index = nullptr;
			cache = nullptr;
			natm = 0;
			nbas = 0;
			nenv = 0;
			nelec = 0;
			nbf_tot = 0;
			nbf_max = 0;
		}

		~IntegralEngineHandle() {
			FREE(atm);
			FREE(bas);
			FREE(env);
			FREE(shell_to_index);
			FREE(cache);
		}
	};

	template<IntegralEngine ENGINE>
	void  initialize_handle (double* molecule, int molecule_size, double* basis, int basis_size, IntegralEngineHandle<ENGINE>& handle);

	template<IntegralEngine ENGINE>
	void  destroy_handle    (IntegralEngineHandle<ENGINE>& handle);

	template<IntegralEngine ENGINE, LACIS_TYPE Type, typename T>
	void cint1e_hf (const IntegralEngineHandle<ENGINE>& handle, LACIS<Type, T>& Ts, LACIS<Type, T>& S);

	template<IntegralEngine ENGINE, LACIS_TYPE Type, typename T>
	void cint2e_hf_eri (const IntegralEngineHandle<ENGINE>& handle, LACIS<Type, T>& ERI);

	template<IntegralEngine ENGINE, LACIS_TYPE Type>
	double run_rhf (const IntegralEngineHandle<ENGINE>& handle);

}

#endif //NEWSCF_HPP
