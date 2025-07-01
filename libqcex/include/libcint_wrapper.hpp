#ifndef NEWSCF_LIBCINT_WRAPPER_HPP
#define NEWSCF_LIBCINT_WRAPPER_HPP

#ifndef CINT_OVERLAP_CACHE
#define CINT_OVERLAP_CACHE 1024
#endif

#include "matrix.hpp"

namespace newscf::qcex {

    using newscf::libmatrix::Matrix;

    void prepare_libcint_data(
        const double* geom,     // [NumAtoms, Z1, X1, Y1, Z1, ...]
        const double* basis,    // [Z, NumShells, L, nprim, exp, coef, ...]
        int basis_size,
        int** atm_out, int* natm_out,
        int** bas_out, int* nbas_out,
        double** env_out, int* nenv_out
    );

    void destroy_environment(int* atm, int* bas, double* env);

    void create_overlap_matrix(int* atm, int natm, int* bas, int nbas, double* env, int nenv, Matrix<double>** overlap_matrix);

    void create_kinetic_matrix(int* atm, int natm, int* bas, int nbas, double* env, int nenv, Matrix<double>** kin_matrix);

    void create_nuclear_matrix(int* atm, int natm, int* bas, int nbas, double* env, int nenv, Matrix<double>** nuc_matrix);


}

#endif