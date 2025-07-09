#include <cstdlib>
#include <iostream>
#include <cstring>
#include <cmath>
#include <functional>

#include "newscf/newscf.hpp"
#include "newscf/ndtx.hpp"
#include "newscf/common.hpp"

extern "C" {
    #include "cint.h"
    #include "cint_funcs.h"
}

using namespace newscf;
using namespace newscf::ndtx;

#ifndef LIBCINT_CACHE_SIZE
#define LIBCINT_CACHE_SIZE 2048*18
#endif

typedef std::function<void (CINTOpt**, int*, int, int*, int, double*)> CINTGenericOptimizer;
typedef std::function<CACHE_SIZE_T (double*, int*, int*, int*, int, int*, int, double*, CINTOpt*, double*)> CINTGenericIntegral;

// TODO; use this later
void build_uacache(const double* geom, int** cache, int natoms) {
    LOG (DEV_INFO, "Building unique atom cache");
    int* tc = reinterpret_cast<int*>(newscf::malloc(sizeof(int) * natoms));
    newscf::memset(*cache, -1, natoms * sizeof(int));
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
    (*cache) = reinterpret_cast<int*>(newscf::malloc(ns * sizeof(int)));
    memcpy(*cache, tc, ns * sizeof(int));
    free(tc);
}

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
    LOG (DEV_INFO, "Scanning atom geometry and finding basis");

    for (int gi = 0, i = 0;
         gi < 4 * natoms;
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
    *atm_out  = reinterpret_cast<int*>(newscf::malloc(sizeof(int) * *natm_out));
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
    *bas_out = reinterpret_cast<int*>(newscf::malloc(sizeof(int) * *nbas_out));
    *nenv_out = env_spaces_needed + PTR_ENV_START;
    LOG(DEV_INFO, "Allocating env array (bytes): "+std::to_string(*nenv_out * sizeof(double)));
    *env_out = reinterpret_cast<double*>(newscf::malloc(sizeof(double) * *nenv_out));
    newscf::memset(*env_out, 0, sizeof(int) * PTR_ENV_START);

    // Populate atm coords to env
    LOG(DEV_INFO, "Loading geometry information to env.");
    for (int i = PTR_ENV_START, gi = 0; gi < 4 * natoms; i += 3, gi += 4) {
        // gi + 0 has Z-val ; not needed
        // libcint expects geometry in Bohr units
        (*env_out) [i + 0] = geom[gi + 1] * angstrom_to_bohr; // X
        (*env_out) [i + 1] = geom[gi + 2] * angstrom_to_bohr; // Y
        (*env_out) [i + 2] = geom[gi + 3] * angstrom_to_bohr; // Z
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
    LOG (DEV_INFO, "(free) Removing libcint atm workspace");
    free(atm);
    LOG (DEV_INFO, "(free) Removing libcint bas workspace");
    free(bas);
    LOG (DEV_INFO, "(free) Removing libcint env workspace");
    free(env);
}

void cint1e_hf (newscf::IntegralEngineHandle* handle, newscf::ndtx::NDTX<double>& T, newscf::ndtx::NDTX<double>& S) {
    double* buffer_K = new double[handle->nbf_max * handle->nbf_max];
    double* buffer_N = new double[handle->nbf_max * handle->nbf_max];
    double* buffer_S = new double[handle->nbf_max * handle->nbf_max];
    int*    shells = new int   [2];
    int*    dims   = new int   [3] {1, 1, 1};
    T.resizeToMatrix(handle->nbf_tot);
    T.fill(0);
    S.resizeToMatrix(handle->nbf_tot);
    T.fill(0);
    CINTOpt* opt_K = nullptr;
    CINTOpt* opt_N = nullptr;
    CINTOpt* opt_S = nullptr;
    int1e_kin_optimizer (&opt_K, handle->atm, handle->natm, handle->bas, handle->nbas, handle->env);
    int1e_nuc_optimizer (&opt_N, handle->atm, handle->natm, handle->bas, handle->nbas, handle->env);
    int1e_ovlp_optimizer(&opt_S, handle->atm, handle->natm, handle->bas, handle->nbas, handle->env);
    for (int i = 0; i < handle->nbas; i++) {
        for (int j = 0; j <= i; j++) {
            shells[0] = i;
            shells[1] = j;
            dims[0]   = CINTcgto_spheric(i, handle->bas);
            dims[1]   = CINTcgto_spheric(j, handle->bas);
            int ioff  = handle->shell_to_index[i];
            int joff  = handle->shell_to_index[j];
            newscf::memset(buffer_K, 0, sizeof(double) * dims[0] * dims[1]);
            newscf::memset(buffer_N, 0, sizeof(double) * dims[0] * dims[1]);
            newscf::memset(buffer_S, 0, sizeof(double) * dims[0] * dims[1]);
            int ncomp_K = int1e_kin_sph  (buffer_K, dims, shells, handle->atm, handle->natm, handle->bas, handle->nbas, handle->env, opt_K, handle->cache);
            int ncomp_N = int1e_nuc_sph  (buffer_N, dims, shells, handle->atm, handle->natm, handle->bas, handle->nbas, handle->env, opt_N, handle->cache);
            int ncomp_S = int1e_ovlp_sph (buffer_S, dims, shells, handle->atm, handle->natm, handle->bas, handle->nbas, handle->env, opt_S, handle->cache);
            if (ncomp_K != 1 || ncomp_N != 1 || ncomp_S != 1) HANDLE_ERROR ("libcint 1e calculation failed.", 101);
            for (int ii = 0; ii < dims[0]; ii++) {
                for (int jj = 0; jj < dims[1]; jj++) {
                    double Kval = buffer_K[ii * dims[1] + jj];
                    double Nval = buffer_N[ii * dims[1] + jj];
                    double Sval = buffer_S[ii * dims[1] + jj];
                    T.matrixSet(ioff + ii, joff + jj, Kval + Nval);
                    T.matrixSet(joff + jj, ioff + ii, Kval + Nval);
                    S.matrixSet(ioff + ii, joff + jj, Sval);
                    S.matrixSet(joff + jj, ioff + ii, Sval);
                }
            }
        }
    }
    CINTdel_optimizer (&opt_K);
    CINTdel_optimizer (&opt_N);
    CINTdel_optimizer (&opt_S);
    delete[] buffer_K;
    delete[] buffer_N;
    delete[] buffer_S;
    delete[] dims;
    delete[] shells;
}

void build_eri_tensor (IntegralEngineHandle* handle, newscf::ndtx::NDTX<double>& ERI) {
	// TODO: Parallelize via OpenMP

    const int nbf = handle->nbf_tot;
    const int nshells = handle->nbas;
    ERI.resizeToTensor4D(nbf);
    ERI.fill(0);

    CINTOpt* opt = nullptr;
    int2e_optimizer(&opt, handle->atm, handle->natm, handle->bas, handle->nbas, handle->env);
    double* cache = new double[LIBCINT_CACHE_SIZE * 5]; // TODO: check PySCF cache sizes

    const int max_prim = handle->nbf_max * handle->nbf_max * handle->nbf_max * handle->nbf_max;
    double* buf = new double[max_prim];

    for (int i = 0; i < nshells; ++i) {
        int di = CINTcgto_spheric(i, handle->bas);
        int ioff = handle->shell_to_index[i];

        for (int j = 0; j <= i; ++j) {
            int dj = CINTcgto_spheric(j, handle->bas);
            int joff = handle->shell_to_index[j];

            for (int k = 0; k < nshells; ++k) {
                int dk = CINTcgto_spheric(k, handle->bas);
                int koff = handle->shell_to_index[k];

                for (int l = 0; l <= k; ++l) {
                    int dl = CINTcgto_spheric(l, handle->bas);
                    int loff = handle->shell_to_index[l];

                    // Enforce (ij) ≤ (kl) in lex order
                    if ((i * nshells + j) < (k * nshells + l)) continue;

                    int shls[4] = {i, j, k, l};
                    int dims[4] = {di, dj, dk, dl};

                    int ncomp = int2e_sph(buf, dims, shls,
                                          handle->atm, handle->natm,
                                          handle->bas, handle->nbas,
                                          handle->env, opt, cache);
                    if (ncomp != 1) HANDLE_ERROR("int2e_sph failed", 1);

                    for (int ii = 0; ii < di; ++ii) {
                        int mu = ioff + ii;
                        for (int jj = 0; jj < dj; ++jj) {
                            int nu = joff + jj;
                            for (int kk = 0; kk < dk; ++kk) {
                                int lam = koff + kk;
                                for (int ll = 0; ll < dl; ++ll) {
                                    int sig = loff + ll;

                                    int idx = (((ii * dj + jj) * dk + kk) * dl + ll);
                                    double val = buf[idx];

                                    ERI.tensor4DSet(mu, nu, lam, sig, val);
                                    ERI.tensor4DSet(nu, mu, lam, sig, val);
                                    ERI.tensor4DSet(mu, nu, sig, lam, val);
                                    ERI.tensor4DSet(nu, mu, sig, lam, val);
                                    ERI.tensor4DSet(lam, sig, mu, nu, val);
                                    ERI.tensor4DSet(lam, sig, nu, mu, val);
                                    ERI.tensor4DSet(sig, lam, mu, nu, val);
                                    ERI.tensor4DSet(sig, lam, nu, mu, val);

                                }
                            }
                        }
                    }
                }
            }
        }
    }

    CINTdel_2e_optimizer(&opt);
    delete[] cache;
    delete[] buf;
}

namespace newscf {

    void  initialize_handle (double* molecule, int molecule_size, double* basis, int basis_size, IntegralEngineHandle* handle) {
        prepare_libcint_data (molecule, molecule_size, basis, basis_size, &(handle->atm), &(handle->natm), &(handle->bas), &(handle->nbas), &(handle->env), &(handle->nenv));
        /*for(int i = 0; i < handle->natm; i++)
            std::cout << handle->atm[i] << "\n";
        std::cout << "\n";
        for(int i = 0; i < handle->nbas; i++)
            std::cout << handle->bas[i] << "\n";
        std::cout << "\n";
        for(int i = 0; i < handle->nenv; i++)
            std::cout << handle->env[i] << "\n";
        std::cout << "\n";*/
        handle->natm /= ATM_SLOTS;
        handle->nbas /= BAS_SLOTS;
        handle->shell_to_index = reinterpret_cast<int*> (malloc(sizeof(int) * (handle->nbas + 1)));
        handle->shell_to_index[0] = 0;
        handle->nbf_tot = 0;
        handle->nbf_max = 0;
        for (int i = 0; i < handle->nbas; i++) {
            int nbf = CINTcgto_spheric(i, handle->bas);  // correct number of contracted AOs
            handle->shell_to_index[i + 1] = handle->shell_to_index[i] + nbf;
            if (nbf > handle->nbf_max) handle->nbf_max = nbf;
            handle->nbf_tot += nbf;
        }
    	for (int i = 0; i < handle->natm; i++) {
    		int Z = handle->atm[i * ATM_SLOTS + ATOM_OF];
    		std::cout << "Z: " << Z << std::endl;
    		handle->nelec += Z;
    	}
        handle->cache = reinterpret_cast<double*> (malloc(sizeof(double) * LIBCINT_CACHE_SIZE));
    }

    void  destroy_handle   (IntegralEngineHandle* handle) {
        destroy_libcint_environment (handle->atm, handle->bas, handle->env);
        free(handle->shell_to_index);
        free(handle->cache);
    }

	void calculate_eri_tensor                  (IntegralEngineHandle* handle, ndtx::NDTX<double>& ERI) {
	    build_eri_tensor(handle, ERI);
    }

    void calculate_hf_matrices                 (IntegralEngineHandle* handle, NDTX<double>& T, NDTX<double>& S) {
        cint1e_hf (handle, T, S);
    }

}