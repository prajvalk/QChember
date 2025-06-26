#ifndef NEWSCF_LIBCINT_WRAPPER_HPP
#define NEWSCF_LIBCINT_WRAPPER_HPP

namespace newscf::qcex {

    void prepare_libcint_data(
        const double* geom,     // [NumAtoms, Z1, X1, Y1, Z1, ...]
        const double* basis,    // [Z, NumShells, L, nprim, exp, coef, ...]
        int** atm_out, int* natm_out,
        int** bas_out, int* nbas_out,
        double** env_out, int* nenv_out
    );

    void destroy_environment(int* atm, int* bas, double* env);

}

#endif