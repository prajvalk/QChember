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
}


#ifndef LIBCINT_CACHE_SIZE
// 128kB Cache
#define LIBCINT_CACHE_SIZE 1024 * 128
#endif

#include "newscf/newscf.hpp"
#include "newscf/newscf_utils.hpp"

namespace newscf {

    void prepare_libcint_data(
    const double* geom, int natoms,
    const double* basis, int basis_size,
    int** atm_out, int* natm_out,
    int** bas_out, int* nbas_out,
    double** env_out, int* nenv_out
) {
    // Compute Sizes
    int* atom_bas_ptr = new int[natoms]; // index-pointers to basis array to where the basis for that atom are present
    int* atom_shells_cont = new int[natoms]; // how many shells are present for each atom
     // env array should hold the atom coords (3 * N_atm; 4-D coordinate atoms are not supported yet)
     //                and basis contraction (expo, coeff); needs to be calculated
    int env_spaces_needed = 3 * natoms;
    DEVLOG (DEV_INFO, "Scanning atom geometry and finding basis");

    for (int gi = 0, i = 0;
         gi < 4 * natoms;
         gi += 4, i++) // for every i-th atom, geom[gi] stores Z-val for Atom, (x, y, z) in (i+1, i+2, i+3)
    {
        // To future self: TODO implement some kind of an atom cache
        // so that if the atom basis has been seen before don't scan basis again
        double atom_z = geom[gi];
        DEVLOG (DEV_INFO, "ATOM "+std::to_string(i));
        for (int j = 0; j < basis_size;) {
            double basis_z = basis[j];
            if (static_cast<int>(atom_z) == static_cast<int>(basis_z)) {
                // We have found a basis for atom in geom
                atom_bas_ptr[i] = j;
                atom_shells_cont[i] = static_cast<int>(basis[j+1]);
                DEVLOG(DEV_INFO, "No. of Shells: "+std::to_string(atom_shells_cont[i]));
                // compute how much space of env we need for this atom
                int shells_to_account = atom_shells_cont[i];
                int k = 0;
                for(k = j + 2; k < basis_size;) {
                    // k-th index is ang. momentum for shell
                    // k+1-th index is contractions
                    int contr = static_cast<int>(basis[k+1]);
                    // we need to spaces in env for each contraction
                    env_spaces_needed += 2 * contr;
                    // account for this shell
                    shells_to_account--;
                    // check if need to account for more
                    if (shells_to_account == 0) break;
                    // skip ahead of contractions
                    k += 2 * contr + 2;
                }
                // we can skip ahead now, since we have read the basis for this atom
                j = k;
                // but then we don't have to look further, now do we?
                break;
            } else {
                // We need to skip ahead of this basis we found
                int shells_to_skip = static_cast<int>(basis[j+1]);
                int k = 0;
                for(k = j+2; k < basis_size;) {
                    // k-th index is angular momentum of shell; we don't need it at this time
                    // k+1-th index is contractions
                    int contr = static_cast<int>(basis[k+1]);
                    // each contraction occupies 2 spaces (expo, coeff); let's skip ahead
                    k += 2 * contr + 2;
                    //std::cout << k << " " << shells_to_skip << "\n";
                    // we should have skipped ahead of this shell
                    shells_to_skip--;
                    // check if we have more shells to skip
                    if (shells_to_skip == 0) break;
                }
                // if everything went right; k will have the next Z-val for next atom
                j = k;
                continue;
            }
            // for future self, TODO: verify if any atoms are missing a basis
        }
    }

    // Allocate Atom Info Array
    *natm_out = ATM_SLOTS * natoms;
    DEVLOG(DEV_INFO, "Scan done. Allocating atm array (bytes): "+std::to_string(*natm_out * sizeof(int)));
    void* atmspace = MALLOC (sizeof(int) * *natm_out);
    *atm_out  = reinterpret_cast<int*>(atmspace);
    int shells_needed = 0;

    // Populate Atom Info
    DEVLOG(DEV_INFO, "Loading atom information to atm.");
    for(int i = 0; i < natoms; i++) {
        // Nuclear Charge; should be the index pointed to in atom_ptr
        (*atm_out) [i * ATM_SLOTS + CHARGE_OF] = static_cast<int>(basis[atom_bas_ptr[i]]);
        // COORD PTR; env will have sequential storage (all atom coords are stored first then basis info is stored)
        (*atm_out) [i * ATM_SLOTS + PTR_COORD] = 3 * i + PTR_ENV_START;
        // Nuclear Model
        (*atm_out) [i * ATM_SLOTS + NUC_MOD_OF] = 1;
        // ECP Zeta pointer (unsupported for now)
        (*atm_out) [i * ATM_SLOTS + PTR_ZETA] = 0;
        // Fraction Charge
        (*atm_out) [i * ATM_SLOTS + PTR_FRAC_CHARGE] = 0;
        // Reserved Slot
        (*atm_out) [i * ATM_SLOTS + RESERVE_ATMSLOT] = 0;
        shells_needed += atom_shells_cont[i];
    }

    // For future self; some memory optimization can be done here, because of redundant atoms
    // Allocate Basis Info and Env Array
    *nbas_out = BAS_SLOTS * shells_needed;
    DEVLOG(DEV_INFO, "Allocating bas array (bytes): "+std::to_string(*nbas_out * sizeof(int)));
    *bas_out = reinterpret_cast<int*>(malloc(sizeof(int) * *nbas_out));
    *nenv_out = env_spaces_needed + PTR_ENV_START;
    DEVLOG(DEV_INFO, "Allocating env array (bytes): "+std::to_string(*nenv_out * sizeof(double)));
    *env_out = reinterpret_cast<double*>(malloc(sizeof(double) * *nenv_out));

    // TODO: Implement Memset
    memset(*env_out, 0, sizeof(int) * PTR_ENV_START);

    // Populate atm coords to env
    DEVLOG(DEV_INFO, "Loading geometry information to env.");
    for (int i = PTR_ENV_START, gi = 0; gi < 4 * natoms; i += 3, gi += 4) {
        // gi + 0 has Z-val ; not needed
        // libcint expects geometry in Bohr units
        (*env_out) [i + 0] = geom[gi + 1] * utils::angstrom_to_bohr; // X
        (*env_out) [i + 1] = geom[gi + 2] * utils::angstrom_to_bohr; // Y
        (*env_out) [i + 2] = geom[gi + 3] * utils::angstrom_to_bohr; // Z
    }

    // Working Implementaion
    int env_ptr = 3 * natoms + PTR_ENV_START;
    int bas_ptr = 0;
    for (int atmi = 0; atmi < natoms; atmi++) {
        int nshells = atom_shells_cont[atmi];
        int basi = atom_bas_ptr[atmi] + 2;
        for (int shelli = 0; shelli < nshells; shelli++) {
            int l = static_cast<int>(basis[basi]);       // Shell ANgular Momentum
            int ncontr = static_cast<int>(basis[basi+1]); // No. of contractions
            basi += 2; // skip ahead of shell header
            int next_shell = basi + 2 * ncontr;
            int env_exp_ptr = env_ptr + 0;
            int env_coef_ptr = env_ptr + ncontr;
            (*bas_out)[bas_ptr + ATOM_OF]   = atmi;       // Atom Index
            (*bas_out)[bas_ptr + ANG_OF]    = l;          // Shell Angular Momentum
            (*bas_out)[bas_ptr + NPRIM_OF]  = ncontr;     // Number of Primitive GTOs
            (*bas_out)[bas_ptr + NCTR_OF]   = 1;          // Compound contractions unsupported
            (*bas_out)[bas_ptr + KAPPA_OF]  = 0;         // Spherical type basis; cartesian or spinor basis will be implemented in the future
            (*bas_out)[bas_ptr + PTR_EXP]   = env_exp_ptr;  // where are my exponents strored in env?
            (*bas_out)[bas_ptr + PTR_COEFF] = env_coef_ptr; // where are my coefficients stored?
            (*bas_out)[bas_ptr + RESERVE_BASLOT] = 0;      // No idea

            // BSE doesn't (always) normalize their CGTOs and PGTOs; and QChem isn't very happy
            // sigh........
            // guess how I found out?
            //
            double* exp = new double[ncontr];  // TODO: Eliminate the need for a seperate workspace
            double* coeff = new double[ncontr];

            // 1. Primitive normalization
            for (int i = 0, j = 0; i < 2 * ncontr; i += 2, j++) {
                exp[j] = basis[basi + i];
                coeff[j] = basis[basi + i + 1] * CINTgto_norm(l, exp[j]);
            }

            // 2. Contracted normalization
            double norm_factor = 0.0;
            for (int i = 0; i < ncontr; i++) {
                for (int j = 0; j < ncontr; j++) {
                    double a = exp[i], b = exp[j];
                    double p = a + b;
                    double Sij = 0.5 * std::tgamma((2 * l + 3) / 2.0) / std::pow(p, (2 * l + 3) / 2.0);
                    norm_factor += coeff[i] * coeff[j] * Sij;
                }
            }
            norm_factor = 1.0 / std::sqrt(norm_factor);

            // 3. Store normalized exponents and coefficients
            for (int i = 0; i < ncontr; i++) {
                (*env_out)[env_exp_ptr++] = exp[i];
                (*env_out)[env_coef_ptr++] = coeff[i] * norm_factor;
            }

            delete[] exp;
            delete[] coeff;
            //
            // BSE; please normalize
            basi = next_shell;
            bas_ptr += BAS_SLOTS;
            env_ptr += 2 * ncontr;
        }
    }

    // Clean Workspace
    delete[] atom_bas_ptr;
    delete[] atom_shells_cont;
}

    void destroy_libcint_environment(int* atm, int* bas, double* env) {
        DEVLOG (DEV_INFO, "(free) Removing libcint atm workspace");
        free(atm);
        DEVLOG (DEV_INFO, "(free) Removing libcint bas workspace");
        free(bas);
        DEVLOG (DEV_INFO, "(free) Removing libcint env workspace");
        free(env);
    }

    template<>
    void  initialize_handle (double* molecule, int molecule_size, double* basis, int basis_size, IntegralEngineHandle<LIBCINT>& handle) {
        prepare_libcint_data (molecule, molecule_size, basis, basis_size, &(handle.atm), &(handle.natm), &(handle.bas), &(handle.nbas), &(handle.env), &(handle.nenv));
        /*for(int i = 0; i < handle->natm; i++)
            std::cout << handle->atm[i] << "\n";
        std::cout << "\n";
        for(int i = 0; i < handle->nbas; i++)
            std::cout << handle->bas[i] << "\n";
        std::cout << "\n";
        for(int i = 0; i < handle->nenv; i++)
            std::cout << handle->env[i] << "\n";
        std::cout << "\n";*/
        handle.natm /= ATM_SLOTS;
        handle.nbas /= BAS_SLOTS;
        handle.shell_to_index = reinterpret_cast<int*> (malloc(sizeof(int) * (handle.nbas + 1)));
        handle.shell_to_index[0] = 0;
        handle.nbf_tot = 0;
        handle.nbf_max = 0;
        for (int i = 0; i < handle.nbas; i++) {
            int nbf = CINTcgto_spheric(i, handle.bas);  // correct number of contracted AOs
            handle.shell_to_index[i + 1] = handle.shell_to_index[i] + nbf;
            if (nbf > handle.nbf_max) handle.nbf_max = nbf;
            handle.nbf_tot += nbf;
        }
        handle.nelec=0;
        for (int i = 0; i < handle.natm; i++) {
            const int Z = handle.atm[i * ATM_SLOTS + ATOM_OF];
            handle.nelec += Z;
        }
        handle.cache = reinterpret_cast<double*> (malloc(sizeof(double) * LIBCINT_CACHE_SIZE));
    }

    template<>
    void  destroy_handle   (IntegralEngineHandle<LIBCINT>& handle) {
    }

}
