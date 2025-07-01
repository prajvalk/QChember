#include <cstdlib>
#include <iostream>
#include <cstring>
#include <cmath>
#include "qcex_utils.hpp"

#include "libcint_wrapper.hpp"
#include "logging_api.hpp"

extern "C" {
    #include "cint.h"
    #include "cint_funcs.h"
}

void _build_uacache(const double* geom, int** cache, int natoms) {
    LOG (DEV_INFO, "Building unique atom cache");
    int* tc = reinterpret_cast<int*>(malloc(sizeof(int) * natoms));
    memset(*cache, -1, natoms * sizeof(int));
    int cptr = 0, ns = 0;
    for(int i = 0; i < natoms; i++) {
        bool unique = true;
        int cval = geom[1 + i * 4];
        for(int j = i - 1; j >= 0; j--)
            if (tc[j] == cval)
                unique = false;
        if(unique) {
            tc[cptr++] = cval;
            ns++;
        }
    }
    (*cache) = reinterpret_cast<int*>(malloc(ns * sizeof(int)));
    memcpy(*cache, tc, ns * sizeof(int));
    free(tc);
}

void newscf::qcex::prepare_libcint_data(
    const double* geom,
    const double* basis, int basis_size,
    int** atm_out, int* natm_out,
    int** bas_out, int* nbas_out,
    double** env_out, int* nenv_out
) {
    int natoms = static_cast<int>(geom[0]);

    // Compute Sizes
    int* atom_bas_ptr = new int[natoms]; // index-pointers to basis array to where the basis for that atom are present
    int* atom_shells_cont = new int[natoms]; // how many shells are present for each atom
     // env array should hold the atom coords (3 * N_atm; 4-D coordinate atoms are not supported yet) 
     //                and basis contraction (expo, coeff); needs to be calculated
    int env_spaces_needed = 3 * natoms;
    LOG (DEV_INFO, "Scanning atom geometry and finding basis");

    for (int gi = 1, i = 0;
         gi < 4 * natoms + 1; 
         gi += 4, i++) // for every i-th atom, geom[gi] stores Z-val for Atom, (x, y, z) in (i+1, i+2, i+3)
    {
        // To future self: TODO implement some kind of an atom cache 
        // so that if the atom basis has been seen before don't scan basis again
        double atom_z = geom[gi];
        LOG (DEV_INFO, "ATOM "+std::to_string(i));
        for (int j = 0; j < basis_size;) {
            double basis_z = basis[j];
            if (static_cast<int>(atom_z) == static_cast<int>(basis_z)) {
                // We have found a basis for atom in geom   
                atom_bas_ptr[i] = j;
                atom_shells_cont[i] = static_cast<int>(basis[j+1]);
                LOG(DEV_INFO, "No. of Shells: "+std::to_string(atom_shells_cont[i]));
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
    LOG(DEV_INFO, "Scan done. Allocating atm array (bytes): "+std::to_string(*natm_out * sizeof(int)));
    *atm_out  = reinterpret_cast<int*>(malloc(sizeof(int) * *natm_out));
    int shells_needed = 0;

    // Populate Atom Info
    LOG(DEV_INFO, "Loading atom information to atm.");
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
    LOG(DEV_INFO, "Allocating bas array (bytes): "+std::to_string(*nbas_out * sizeof(int)));
    *bas_out = reinterpret_cast<int*>(malloc(sizeof(int) * *nbas_out));
    *nenv_out = env_spaces_needed + PTR_ENV_START;
    LOG(DEV_INFO, "Allocating env array (bytes): "+std::to_string(*nenv_out * sizeof(double)));
    *env_out = reinterpret_cast<double*>(malloc(sizeof(double) * *nenv_out));
    memset(*env_out, 0, sizeof(int) * PTR_ENV_START);

    // Populate atm coords to env
    LOG(DEV_INFO, "Loading geometry information to env.");
    for (int i = PTR_ENV_START, gi = 1; gi < 4 * natoms + 1; i += 3, gi += 4) {
        // gi + 0 has Z-val ; not needed
        // libcint expects geometry in Bohr units
        (*env_out) [i + 0] = geom[gi + 1] * angstrom_to_bohr; // X
        (*env_out) [i + 1] = geom[gi + 2] * angstrom_to_bohr; // Y
        (*env_out) [i + 2] = geom[gi + 3] * angstrom_to_bohr; // Z
    }

    // More efficient to populate env and bas simultaneously for basis info
    // Obselete & Buggy Routine; see below for current implementation
    /*int env_offset = 3 * natoms + PTR_ENV_START; // skip coordinate info in env
    int env_ptr = env_offset;
    int bas_ptr = 0;
    LOG(DEV_INFO, "Loading basis information to bas and env.");
    for (int atm_i = 0; atm_i < natoms; atm_i++) {
        LOG(DEV_INFO, "Loading basis information for ATOM "+std::to_string(atm_i));
        int atm_shell_no = atom_shells_cont[atm_i]; // get number of shells
        int atm_bas_i = atom_bas_ptr[atm_i] + 2; // atom_bas_ptr at [atm_i] and [atm_i+1] will have Z-val and Num of Shells; not needed
        for (int shell_i = 0; shell_i < atm_shell_no; shell_i++) {
            // Load contractions first
            int shell_l = static_cast<int>(basis[atm_bas_i]); // angular momentum info
            int shell_contr = static_cast<int>(basis[atm_bas_i+1]); // number of contractions for this shell
            for (int j = atm_bas_i+2, i = 0; j < atm_bas_i+2*shell_contr; j += 2, i += 2) {
                (*env_out)[env_ptr + i + 0] = basis[j];              // exponents in contigous memory segment
                (*env_out)[env_ptr + i + shell_contr] = basis[j+1];  // coeffecients in contigous memory segement offset by number of contractions
            }
            (*bas_out)[bas_ptr + ATOM_OF]   = atm_i;       // Atom Index
            (*bas_out)[bas_ptr + ANG_OF]    = shell_l;     // Shell Angular Momentum
            (*bas_out)[bas_ptr + NPRIM_OF]  = shell_contr; // Number of Primitive GTOs
            (*bas_out)[bas_ptr + NCTR_OF]   = 1;           // Compound contractions unsupported
            (*bas_out)[bas_ptr + KAPPA_OF]  = -1;          // Spherical type basis; cartesian or spinor basis will be implemented in the future
            (*bas_out)[bas_ptr + PTR_EXP]   = env_ptr + 0; // where are my exponents strored in env?
            (*bas_out)[bas_ptr + PTR_COEFF] = env_ptr + shell_contr; // where are my coefficients stored?
            (*bas_out)[bas_ptr + RESERVE_BASLOT] = 0;      // No idea
            bas_ptr += BAS_SLOTS;
            env_ptr += 2 * shell_contr;
            atm_bas_i += 2 * shell_contr + 2; // move to the next shell
        }
    }*/

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

    // (In)sanity Check
    /*if (env_ptr >= *nenv_out) {
        // lmao, ded :(
        // it's pointless anyway; would segfault instantly much before
    }*/

    // Checking routines
    // Atom Check
    /*std::cout << "Atom Geometry Check \n";
    for (int i = 0; i < natoms; i++) {
        int zval = (*atm_out)[i * ATM_SLOTS + CHARGE_OF];
        int coordptr = (*atm_out)[i * ATM_SLOTS + PTR_COORD];
        std::cout << "Z:" << get_atomic_symbol(zval) << " X:" << (*env_out)[coordptr] << " Y:" << (*env_out)[coordptr + 1] << " Z:" << (*env_out)[coordptr + 2] << "\n";
    }

    // Basis Check
    std::cout << "Basis Check \n";
    for (int i = 0; i < *nbas_out; i += BAS_SLOTS) {
        int atmi = (*bas_out)[i + ATOM_OF];
        int atmz = (*atm_out)[atmi * ATM_SLOTS];
        int shell_l = (*bas_out)[i + ANG_OF];
        int ncontr = (*bas_out)[i + NPRIM_OF];
        int env_exp_ptr = (*bas_out)[i + PTR_EXP];
        int env_coeff_ptr = (*bas_out)[i + PTR_COEFF];
        std::cout << get_atomic_symbol(atmz) << "/" << shell_l << "\n";
        std::cout << "Exponents: \n";
        for (int j = 0; j < ncontr; j++)
            std::cout << (*env_out)[env_exp_ptr + j] << "\n";
        std::cout << "\n";
        std::cout << "Coefficients: \n";
        for (int j = 0; j < ncontr; j++)
            std::cout << (*env_out)[env_coeff_ptr + j] << "\n";
        std::cout << "\n";
    }*/

    // Clean Workspace
    delete[] atom_bas_ptr;
    delete[] atom_shells_cont;
}


void newscf::qcex::destroy_environment(int* atm, int* bas, double* env) {
    LOG (DEV_INFO, "(free) Removing libcint atm workspace");
    free(atm);
    LOG (DEV_INFO, "(free) Removing libcint bas workspace");
    free(bas);
    LOG (DEV_INFO, "(free) Removing libcint env workspace");
    free(env);
}

void newscf::qcex::create_overlap_matrix(
    int* atm, int natm,
    int* bas, int nbas,
    double* env, int nenv,
    Matrix<double>** overlap_matrix  // caller must pass address of pointer
) {
    int natoms = natm / ATM_SLOTS;
    int nshells = nbas / BAS_SLOTS;
    int nbf_tot = 0;

    double* cache = new double[CINT_OVERLAP_CACHE];
    CINTOpt* opt = nullptr;
    int1e_ovlp_optimizer(&opt, atm, natm, bas, nbas, env);

    int* shell_to_index = new int[nshells + 1];
    shell_to_index[0] = 0; 
    int buffmax = 0;
    for (int i = 0; i < nshells; i++) {
        int l = static_cast<int>(bas[BAS_SLOTS * i + ANG_OF]);
        int nctr = static_cast<int>(bas[BAS_SLOTS * i + NPRIM_OF]);
        int nbf = (2 * l + 1) * nctr;
        shell_to_index[i + 1] = shell_to_index[i] + nbf;
        nbf_tot += nbf;
        if (buffmax < nbf * nbf) buffmax = nbf * nbf;
    }

    (*overlap_matrix) = new Matrix<double>(nbf_tot, nbf_tot);

    double* buf = reinterpret_cast<double*>(malloc(sizeof(double) * buffmax));
    int* shells = new int[2];
    int* dims = new int[3]{1, 1, 1};

    for (int i = 0; i < nshells; i++) {
        for (int j = 0; j <= i; j++) {
            shells[0] = i;
            shells[1] = j;
            int di = CINTcgto_spheric(i, bas);
            int dj = CINTcgto_spheric(j, bas);
            dims[0] = di;
            dims[1] = dj;

            memset(buf, 0, sizeof(double) * di * dj);

            int ncomp = int1e_ovlp_sph(buf, dims, shells, atm, natm, bas, nbas, env, opt, cache);

            if (ncomp != 1) {
                LOG(DEV_ERROR, "libcint return value is zero for overlap matrix computation.");
                HANDLE_ERROR ("Internal error in computing the overlap matrix.", 561);
            }

            int ioff = shell_to_index[i];
            int joff = shell_to_index[j];

            for (int ii = 0; ii < di; ii++) {
                for (int jj = 0; jj < dj; jj++) {
                    double val = buf[ii * dj + jj];
                    (*overlap_matrix)->set(ioff + ii, joff + jj, val);
                    (*overlap_matrix)->set(joff + jj, ioff + ii, val);
                }
            }
        }
    }

    CINTdel_optimizer(&opt);

    delete[] shell_to_index;
    delete[] cache;
    delete[] shells;
    delete[] dims;

    free(buf);
}

void newscf::qcex::create_kinetic_matrix(
    int* atm, int natm,
    int* bas, int nbas,
    double* env, int nenv,
    Matrix<double>** out  // caller must pass address of pointer
) {
    int natoms = natm / ATM_SLOTS;
    int nshells = nbas / BAS_SLOTS;
    int nbf_tot = 0;

    double* cache = new double[CINT_OVERLAP_CACHE];
    CINTOpt* opt = nullptr;
    int1e_kin_optimizer(&opt, atm, natoms, bas, nshells, env);

    int* shell_to_index = new int[nshells + 1];
    shell_to_index[0] = 0; 
    int buffmax = 0;
    for (int i = 0; i < nshells; i++) {
        int l = static_cast<int>(bas[BAS_SLOTS * i + ANG_OF]);
        int nctr = static_cast<int>(bas[BAS_SLOTS * i + NPRIM_OF]);
        int nbf = (2 * l + 1) * nctr;
        shell_to_index[i + 1] = shell_to_index[i] + nbf;
        nbf_tot += nbf;
        if (buffmax < nbf * nbf) buffmax = nbf * nbf;
    }

    (*out) = new Matrix<double>(nbf_tot, nbf_tot);

    double* buf = reinterpret_cast<double*>(malloc(sizeof(double) * buffmax));
    int* shells = new int[2];
    int* dims = new int[3]{1, 1, 1};

    for (int i = 0; i < nshells; i++) {
        for (int j = 0; j <= i; j++) {
            shells[0] = i;
            shells[1] = j;
            int di = CINTcgto_spheric(i, bas);
            int dj = CINTcgto_spheric(j, bas);
            dims[0] = di;
            dims[1] = dj;

            memset(buf, 0, sizeof(double) * di * dj);

            int ncomp = int1e_kin_sph(buf, dims, shells, atm, natoms, bas, nshells, env, opt, cache);

            if (ncomp != 1) {
                LOG(DEV_ERROR, "libcint return value is zero for kineric matrix computation.");
                HANDLE_ERROR ("Internal error in computing the kinetic matrix.", 561);
            }

            int ioff = shell_to_index[i];
            int joff = shell_to_index[j];

            for (int ii = 0; ii < di; ii++) {
                for (int jj = 0; jj < dj; jj++) {
                    double val = buf[ii * dj + jj];
                    (*out)->set(ioff + ii, joff + jj, val);
                    (*out)->set(joff + jj, ioff + ii, val);
                }
            }
        }
    }

    CINTdel_optimizer(&opt);

    delete[] shell_to_index;
    delete[] cache;
    delete[] shells;
    delete[] dims;

    free(buf);
}

void newscf::qcex::create_nuclear_matrix(
    int* atm, int natm,
    int* bas, int nbas,
    double* env, int nenv,
    Matrix<double>** out  // caller must pass address of pointer
) {
    int natoms = natm / ATM_SLOTS;
    int nshells = nbas / BAS_SLOTS;
    int nbf_tot = 0;

    double* cache = new double[CINT_OVERLAP_CACHE];
    CINTOpt* opt = nullptr;
    int1e_nuc_optimizer(&opt, atm, natoms, bas, nshells, env);

    int* shell_to_index = new int[nshells + 1];
    shell_to_index[0] = 0; 
    int buffmax = 0;
    for (int i = 0; i < nshells; i++) {
        int l = static_cast<int>(bas[BAS_SLOTS * i + ANG_OF]);
        int nctr = static_cast<int>(bas[BAS_SLOTS * i + NPRIM_OF]);
        int nbf = (2 * l + 1) * nctr;
        shell_to_index[i + 1] = shell_to_index[i] + nbf;
        nbf_tot += nbf;
        if (buffmax < nbf * nbf) buffmax = nbf * nbf;
    }

    (*out) = new Matrix<double>(nbf_tot, nbf_tot);

    double* buf = reinterpret_cast<double*>(malloc(sizeof(double) * buffmax));
    int* shells = new int[2];
    int* dims = new int[3]{1, 1, 1};

    for (int i = 0; i < nshells; i++) {
        for (int j = 0; j <= i; j++) {
            shells[0] = i;
            shells[1] = j;
            int di = CINTcgto_spheric(i, bas);
            int dj = CINTcgto_spheric(j, bas);
            dims[0] = di;
            dims[1] = dj;

            memset(buf, 0, sizeof(double) * di * dj);

            int ncomp = int1e_nuc_sph(buf, dims, shells, atm, natoms, bas, nshells, env, opt, cache);

            if (ncomp != 1) {
                LOG(DEV_ERROR, "libcint return value is zero for nuclear attraction matrix computation.");
                HANDLE_ERROR ("Internal error in computing the nuclear attraction matrix.", 561);
            }

            int ioff = shell_to_index[i];
            int joff = shell_to_index[j];

            for (int ii = 0; ii < di; ii++) {
                for (int jj = 0; jj < dj; jj++) {
                    double val = buf[ii * dj + jj];
                    (*out)->set(ioff + ii, joff + jj, val);
                    (*out)->set(joff + jj, ioff + ii, val);
                }
            }
        }
    }

    CINTdel_optimizer(&opt);

    delete[] shell_to_index;
    delete[] cache;
    delete[] shells;
    delete[] dims;

    free(buf);
}