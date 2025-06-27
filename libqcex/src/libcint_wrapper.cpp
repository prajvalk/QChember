#include <cstdlib>
#include <iostream>
#include <cstring>
#include "qcex_utils.hpp"

#include "libcint_wrapper.hpp"
#include "logging_api.hpp"

extern "C" {
    #include "cint.h"
    #include "cint_funcs.h"
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
    for (int i = 1;
         i < 4 * natoms + 1; 
         i += 4) // for every i-th atom, geom[i] stores Z-val for Atom, (x, y, z) in (i+1, i+2, i+3)
    {
        // To future self: TODO implement some kind of an atom cache 
        // so that if the atom basis has been seen before don't scan basis again
        double atom_z = geom[i];
        for (int j = 0; j < basis_size; j++) {
            double basis_z = basis[j];
            if (static_cast<int>(atom_z) == static_cast<int>(basis_z)) {
                // We have found a basis for atom in geom   
                atom_bas_ptr[i] = j;
                atom_shells_cont[i] = static_cast<int>(basis[j+1]);
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
                    k += 2 * contr;
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
    *atm_out  = reinterpret_cast<int*>(malloc(sizeof(int) * *natm_out));
    int shells_needed = 0;

    // Populate Atom Info
    for(int i = 0; i < *natm_out; i++) {
        // Nuclear Charge; should be the index pointed to in atom_ptr
        *atm_out [i * ATM_SLOTS + CHARGE_OF] = static_cast<int>(basis[atom_bas_ptr[i]]);
        // COORD PTR; env will have sequential storage (all atom coords are stored first then basis info is stored)
        *atm_out [i * ATM_SLOTS + PTR_COORD] = 3 * i;
        // Nuclear Model
        *atm_out [i * ATM_SLOTS + NUC_MOD_OF] = 0;
        // ECP Zeta pointer (unsupported for now)
        *atm_out [i * ATM_SLOTS + PTR_ZETA] = 0;
        // Fraction Charge
        *atm_out [i * ATM_SLOTS + FRAC_CHARGE_NUC] = 0;
        // Reserved Slot
        *atm_out [i * ATM_SLOTS + RESERVE_ATMSLOT] = 0;
        shells_needed += atom_shells_cont[i];
    }

    // For future self; some memory optimization can be done here, because of redundant atoms
    // Allocate Basis Info and Env Array
    *nbas_out = BAS_SLOTS * shells_needed;
    *bas_out = reinterpret_cast<int*>(malloc(sizeof(int) * *nbas_out));
    *nenv_out = env_spaces_needed;
    *env_out = reinterpret_cast<double*>(malloc(sizeof(double) * *nenv_out));

    // Populate atm coords to env
    for (int i = 0, gi = 1; i < 3 * natoms, gi < 4 * natoms + 1; i += 3, gi += 4) {
        // gi + 0 has Z-val ; not needed
        *env_out [i + 0] = geom[gi + 1]; // X
        *env_out [i + 1] = geom[gi + 2]; // Y
        *env_out [i + 1] = geom[gi + 3]; // Z
    }


    // More efficient to populate env and bas simultaneously for basis info
    int env_offset = 3 * natoms; // skip coordinate info in env
    int env_ptr = env_offset;
    int bas_ptr = 0;
    for (int atm_i = 0; atm_i < natoms; atm_i++) {
        int atm_shell_no = atom_shells_cont[atm_i]; // get number of shells
        int atm_bas_i = atom_bas_ptr[atm_i] + 2; // atom_bas_ptr at [atm_i] and [atm_i+1] will have Z-val and Num of Shells; not needed
        for (int shell_i = 0; shell_i < atm_shell_no; shell_i) {
            // Load contractions first
            int shell_l = static_cast<int>(basis[atm_bas_i]); // angular momentum info
            int shell_contr = static_cast<int>(basis[atm_bas_i+1]); // number of contractions for this shell
            for (int j = atm_bas_i+2, i = 0; j < atm_bas_i+2*shell_contr; j += 2, i += 2) {
                *env_out[env_ptr + i + 0] = basis[j];              // exponents in contigous memory segment
                *env_out[env_ptr + i + shell_contr] = basis[j+1];  // coeffecients in contigous memory segement offset by number of contractions
            }
            *bas_out[bas_ptr + ATOM_OF]   = atm_i;       // Atom Index
            *bas_out[bas_ptr + ANG_OF]    = shell_l;     // Shell Angular Momentum
            *bas_out[bas_ptr + NPRIM_OF]  = shell_contr; // Number of Primitive GTOs
            *bas_out[bas_ptr + NCTR_OF]   = 1;           // Compound contractions unsupported
            *bas_out[bas_ptr + KAPPA_OF]  = -1;          // Spherical type basis; cartesian or spinor basis will be implemented in the future
            *bas_out[bas_ptr + PTR_EXP]   = env_ptr + 0; // where are my exponents strored in env?
            *bas_out[bas_ptr + PTR_COEFF] = env_ptr + shell_contr; // where are my coefficients stored?
            *bas_out[bas_ptr + RESERVE_BASLOT] = 0;      // No idea
            bas_ptr += BAS_SLOTS;
            env_ptr += 2 * shell_contr;
        }
    }

    // (In)sanity Check
    if (env_ptr >= *nenv_out) {
        // lmao, ded :(
        // it's pointless anyway; would segfault instantly much before
    }

    // Clean Workspace
    delete[] atom_bas_ptr;
    delete[] atom_shells_cont;

    /*const double* gptr = &geom[1];

    // === FIRST PASS: count shells and env size ===
    int total_shells = 0;
    int total_env = 3 * natoms; // space for atomic coordinates
    const double* bptr = basis;

    while (*bptr != 0.0) {
        int Z = static_cast<int>(bptr[0]);
        int nshell = static_cast<int>(bptr[1]);
        bptr += 2;

        for (int i = 0; i < nshell; ++i) {
            int L     = static_cast<int>(bptr[0]);
            int nprim = static_cast<int>(bptr[1]);
            bptr += 2 + 2 * nprim;

            total_shells++;
            total_env += 2 * nprim; // exponents and coefficients
        }
    }

    *natm_out = natoms;
    *nbas_out = total_shells;
    *nenv_out = total_env;

    // === Allocate memory ===
    *atm_out = (int*)malloc(sizeof(int) * ATM_SLOTS * natoms);
    *bas_out = (int*)malloc(sizeof(int) * BAS_SLOTS * total_shells);
    *env_out = (double*)malloc(sizeof(double) * total_env);

    int* atm = *atm_out;
    int* bas = *bas_out;
    double* env = *env_out;

    int env_pos = 0;
    int shell_index = 0;

    // === Load coordinates into atm/env ===
    for (int a = 0; a < natoms; ++a) {
        int Z = static_cast<int>(gptr[0]);
        double x = gptr[1], y = gptr[2], z = gptr[3];

        int offset = a * ATM_SLOTS;
        atm[offset + CHARGE_OF]       = Z;
        atm[offset + PTR_COORD]       = env_pos;
        atm[offset + NUC_MOD_OF]      = 0;
        atm[offset + PTR_EXP]         = -1;
        atm[offset + PTR_FRAC_CHARGE] = 0;

        env[env_pos++] = x;
        env[env_pos++] = y;
        env[env_pos++] = z;

        gptr += 4;
    }

    // === Load basis functions for each atom ===
    gptr = &geom[1];  // reset pointer

    for (int a = 0; a < natoms; ++a) {
        int atom_Z = static_cast<int>(gptr[0]);

        // Find corresponding basis block
        const double* bscan = basis;
        bool found = false;

        while (*bscan != 0.0) {
            int bZ = static_cast<int>(bscan[0]);
            int nshell = static_cast<int>(bscan[1]);
            bscan += 2;

            if (bZ != atom_Z) {
                // skip this block
                for (int i = 0; i < nshell; ++i) {
                    int nprim = static_cast<int>(bscan[1]);
                    bscan += 2 + 2 * nprim;
                }
                continue;
            }

            found = true;

            for (int s = 0; s < nshell; ++s) {
                int L     = static_cast<int>(bscan[0]);
                int nprim = static_cast<int>(bscan[1]);
                bscan += 2;

                int ptr_exp  = env_pos;
                for (int i = 0; i < nprim; ++i)
                    env[env_pos++] = bscan[2 * i];
                int ptr_coef = env_pos;
                for (int i = 0; i < nprim; ++i)
                    env[env_pos++] = bscan[2 * i + 1];
                bscan += 2 * nprim;

                int* b = &bas[shell_index * BAS_SLOTS];
                b[ATOM_OF]     = a;
                b[ANG_OF]      = L;
                b[NPRIM_OF]    = nprim;
                b[NCTR_OF]     = 1;
                b[KAPPA_OF]    = -1;
                b[PTR_EXP]     = ptr_exp;
                b[PTR_COEFF]   = ptr_coef;
                b[RESERVE_BASLOT] = 0;

                shell_index++;
            }

            break; // only load one matching basis block
        }

        if (!found) {
            std::cerr << "ERROR: No basis found for atom Z=" << atom_Z << "\n";
            std::exit(1003);
        }

        gptr += 4;
    }

    if (env_pos > total_env) {
        std::cerr << "BUG: env buffer overflow (used=" << env_pos
                  << ", allocated=" << total_env << ")\n";
        std::exit(1004);
    }

    *nenv_out = env_pos; // set exact number used*/
}


void newscf::qcex::destroy_environment(int* atm, int* bas, double* env) {
    LOG (DEV_INFO, "(free) Removing libcint atm workspace");
    free(atm);
    LOG (DEV_INFO, "(free) Removing libcint bas workspace");
    free(bas);
    LOG (DEV_INFO, "(free) Removing libcint env workspace");
    free(env);
}

#include <cstring> // for memset
#include <cassert>
#include <vector>

#define CINT_CACHE_SIZE 10240
void newscf::qcex::create_overlap_matrix(
    int* atm, int natm,
    int* bas, int nbas,
    double* env, int nenv,
    Matrix<double>* overlap_matrix  // caller must pass address of pointer
) {
    // === Sanity Checks ===
    if (!atm || !bas || !env) {
        HANDLE_ERROR("Null pointer in atm, bas, or env", 1001);
    }
    if (natm <= 0 || nbas <= 0 || nenv <= 0) {
        HANDLE_ERROR("Invalid array dimensions in libcint input", 1002);
    }

    LOG(DEV_INFO, "Validating libcint input arrays...");
    for (int a = 0; a < natm; ++a) {
        int charge = atm[a * ATM_SLOTS + CHARGE_OF];
        int coord_ptr = atm[a * ATM_SLOTS + PTR_COORD];
        if (coord_ptr < 0 || coord_ptr + 2 >= nenv) {
            HANDLE_ERROR("Invalid coordinate pointer for atom " + std::to_string(a), 1003);
        }
        if (charge <= 0) {
            HANDLE_ERROR("Invalid atomic number at atom " + std::to_string(a), 1004);
        }
    }

    for (int b = 0; b < nbas; ++b) {
        int L = bas[b * BAS_SLOTS + ANG_OF];
        int nprim = bas[b * BAS_SLOTS + NPRIM_OF];
        int nctr = bas[b * BAS_SLOTS + NCTR_OF];
        int exp_ptr = bas[b * BAS_SLOTS + PTR_EXP];
        int coef_ptr = bas[b * BAS_SLOTS + PTR_COEFF];

        if (L < 0 || L > 5) {
            HANDLE_ERROR("Unsupported angular momentum L=" + std::to_string(L), 1005);
        }
        if (nprim <= 0 || nctr <= 0) {
            HANDLE_ERROR("Invalid nprim or nctr in shell " + std::to_string(b), 1006);
        }
        if (exp_ptr < 0 || coef_ptr < 0 || coef_ptr + nprim > nenv || exp_ptr + nprim > nenv) {
            HANDLE_ERROR("Invalid env pointers in shell " + std::to_string(b), 1007);
        }
    }

    // === Precompute shell->basis offset ===
    int* shell_to_bas = new int[nbas + 1];
    shell_to_bas[0] = 0;
    for (int s = 1; s <= nbas; ++s) {
        int L_prev = bas[(s - 1) * BAS_SLOTS + ANG_OF];
        int nctr_prev = bas[(s - 1) * BAS_SLOTS + NCTR_OF];
        shell_to_bas[s] = shell_to_bas[s - 1] + (2 * L_prev + 1) * nctr_prev;
    }

    int nbasis = shell_to_bas[nbas];
    overlap_matrix = new Matrix<double>(nbasis, nbasis);

    // === Prepare libcint optimizer ===
    CINTOpt* opt = nullptr;
    int1e_ovlp_optimizer(&opt, atm, natm, bas, nbas, env);

    int dims[3] = {0};        // Not used in 1e integrals
    double cache[1024] = {0};
    int shls[2];

    // === Loop over shell pairs ===
    for (int i = 0; i < nbas; ++i) {
        for (int j = 0; j <= i; ++j) {
            shls[0] = i;
            shls[1] = j;

            int Li = bas[i * BAS_SLOTS + ANG_OF];
            int Lj = bas[j * BAS_SLOTS + ANG_OF];
            int ni = bas[i * BAS_SLOTS + NCTR_OF];
            int nj = bas[j * BAS_SLOTS + NCTR_OF];

            int di = (2 * Li + 1) * ni;
            int dj = (2 * Lj + 1) * nj;

            double* buf = new double[di * dj];

            int ncomputed = int1e_ovlp_sph(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache);

            if (ncomputed != di * dj) {
                delete[] buf;
                delete[] shell_to_bas;
                CINTdel_optimizer(&opt);
                HANDLE_ERROR("Integral dimension mismatch: expected " + std::to_string(di * dj) +
                             ", got " + std::to_string(ncomputed), 1008);
            }

            int offset_i = shell_to_bas[i];
            int offset_j = shell_to_bas[j];

            for (int ii = 0; ii < di; ++ii) {
                for (int jj = 0; jj < dj; ++jj) {
                    double val = buf[ii * dj + jj];
                    (overlap_matrix)->set(offset_i + ii, offset_j + jj, val);
                    if (i != j) {
                        (overlap_matrix)->set(offset_j + jj, offset_i + ii, val);
                    }
                }
            }

            delete[] buf;
        }
    }

    delete[] shell_to_bas;
    CINTdel_optimizer(&opt);
    LOG(DEV_INFO, "Overlap matrix computed successfully with size: " +
                  std::to_string(nbasis) + "x" + std::to_string(nbasis));
}
