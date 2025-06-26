#include <cstdlib>
#include <iostream>

#include "libcint_wrapper.hpp"
#include "logging_api.hpp"

#define ATM_SLOTS 6
#define BAS_SLOTS 8

#define ATM_CHARGE 0
#define ATM_COORD 1
#define ATM_NUC_MOD 2
#define ATM_RESERVED 3
#define ATM_POINTER 4
#define ATM_INDEX 5

#define BAS_NPRIM 0
#define BAS_NCTR 1
#define BAS_KAPPA 2
#define BAS_L 3
#define BAS_PTR_EXP 4
#define BAS_PTR_COEFF 5
#define BAS_ATOMIDX 6
#define BAS_ENVIDX 7

void newscf::qcex::prepare_libcint_data(
    const double* geom,
    const double* basis,
    int** atm_out, int* natm_out,
    int** bas_out, int* nbas_out,
    double** env_out, int* nenv_out)
{
    int natoms = static_cast<int>(geom[0]);
    const double* gptr = &geom[1];
    const double* bptr = basis;

    // === FIRST PASS: validate and count ===
    int total_shells = 0;
    int total_env = 0;
    int total_basis_entries = 0;

    // For each atom in geom, ensure basis exists
    for (int a = 0; a < natoms; ++a) {
        int Z = static_cast<int>(gptr[0]);
        const double* bscan = basis;
        bool found = false;

        while (*bscan) {
            int bZ = static_cast<int>(bscan[0]);
            int nshell = static_cast<int>(bscan[1]);
            bscan += 2;
            for (int i = 0; i < nshell; ++i) {
                int nprim = static_cast<int>(bscan[1]);
                bscan += 2 + 2 * nprim;
            }

            if (bZ == Z) {
                found = true;
                break;
            }
        }

        if (!found) {
            std::cerr << "Basis for atomic number Z=" << Z << " not found.\n";
            std::exit(1001);
        }

        gptr += 4;
    }

    // Second pass: compute total shells and env size
    const double* bscan = basis;
    while (*bscan) {
        int nshell = static_cast<int>(bscan[1]);
        bscan += 2;
        for (int i = 0; i < nshell; ++i) {
            int nprim = static_cast<int>(bscan[1]);
            total_shells++;
            total_env += 2 * nprim;
            bscan += 2 + 2 * nprim;
        }
    }

    total_env += 3 * natoms; // xyz coordinates
    *natm_out = natoms;
    *nbas_out = total_shells;
    *nenv_out = total_env;

    long atmmembytes = sizeof(int) * natoms * ATM_SLOTS;
    long basmembytes = sizeof(int) * total_shells * BAS_SLOTS;
    long envmembytes = sizeof(double) * total_env;

    LOG (DEV_INFO, "(malloc) Allocating memory for libcint atm workspace (bytes): "+std::to_string(atmmembytes));
    *atm_out = reinterpret_cast<int*>    (malloc(atmmembytes));
    LOG (DEV_INFO, "(malloc) Allocating memory for libcint bas workspace (bytes): "+std::to_string(basmembytes));
    *bas_out = reinterpret_cast<int*>    (malloc(basmembytes));
    LOG (DEV_INFO, "(malloc) Allocating memory for libcint env workspace (bytes): "+std::to_string(envmembytes));
    *env_out = reinterpret_cast<double*> (malloc(envmembytes));

    *atm_out = new int[natoms * ATM_SLOTS]();
    *bas_out = new int[total_shells * BAS_SLOTS]();
    *env_out = new double[total_env]();

    int* atm = *atm_out;
    int* bas = *bas_out;
    double* env = *env_out;

    int epos = 0;
    gptr = &geom[1];
    int shell_index = 0;

    for (int a = 0; a < natoms; ++a) {
        int Z = static_cast<int>(gptr[0]);
        double x = gptr[1], y = gptr[2], z = gptr[3];
        int atm_offset = a * ATM_SLOTS;

        atm[atm_offset + ATM_CHARGE] = Z;
        atm[atm_offset + ATM_COORD] = epos;
        env[epos++] = x;
        env[epos++] = y;
        env[epos++] = z;

        // Load basis data for this atom
        const double* bscan = basis;
        while (*bscan) {
            int bZ = static_cast<int>(bscan[0]);
            int nshell = static_cast<int>(bscan[1]);
            bscan += 2;

            if (bZ != Z) {
                for (int i = 0; i < nshell; ++i) {
                    int nprim = static_cast<int>(bscan[1]);
                    bscan += 2 + 2 * nprim;
                }
                continue;
            }

            for (int s = 0; s < nshell; ++s) {
                int L = static_cast<int>(bscan[0]);
                int nprim = static_cast<int>(bscan[1]);
                bscan += 2;

                int ptr_exp = epos;
                for (int j = 0; j < nprim; ++j)
                    env[epos++] = bscan[2 * j];
                int ptr_coef = epos;
                for (int j = 0; j < nprim; ++j)
                    env[epos++] = bscan[2 * j + 1];
                bscan += 2 * nprim;

                int* b = &bas[shell_index * BAS_SLOTS];
                b[BAS_NPRIM] = nprim;
                b[BAS_NCTR] = 1;
                b[BAS_KAPPA] = -1;
                b[BAS_L] = L;
                b[BAS_PTR_EXP] = ptr_exp;
                b[BAS_PTR_COEFF] = ptr_coef;
                b[BAS_ATOMIDX] = a;
                b[BAS_ENVIDX] = 0;

                ++shell_index;
            }

            break; // only one basis match per atom
        }

        gptr += 4;
    }

    if (epos > total_env) {
        std::cerr << "BUG: env overflow\n";
        std::exit(1002);
    }

    *nenv_out = epos; // exact number used
}

void newscf::qcex::destroy_environment(int* atm, int* bas, double* env) {
    LOG (DEV_INFO, "(free) Removing libcint atm workspace");
    free(atm);
    LOG (DEV_INFO, "(free) Removing libcint bas workspace");
    free(bas);
    LOG (DEV_INFO, "(free) Removing libcint env workspace");
    free(env);
}