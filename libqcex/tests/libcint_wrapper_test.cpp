#include "libcint_wrapper.hpp"
#include "qcex.hpp"
#include "testing_api.hpp"
#include <cstring>

using namespace newscf::testing;
using namespace newscf::qcex;

int main() {
    TestHandle h;
    h.testsuite = "qcex_libcint_wrapper";

    init_tests();

    add_test("water/def2-TZVPP", h);

    double* geom;
    double* basis;

    load_molecule("test_water.xyz", &geom);
    int basis_size = load_basis("def2-SVP.bas", &basis);

    int* atm;
    int* bas;
    int natm, nbas;
    double* env;
    int nenv;

    prepare_libcint_data (geom, basis, basis_size, &atm, &natm, &bas, &nbas, &env, &nenv);

    Matrix<double> overlap;
    Matrix<double>* matptr = &overlap;
    create_overlap_matrix (atm, natm, bas, nbas, env, nenv, &matptr);

    for (int i = 0; i < matptr->sz_rows; i++) {
        for (int j = 0; j < matptr->sz_cols; j++) {
            if(matptr->get(i, j) != 0) std::cout << i << ", " << j << ", " << matptr->get(i, j) << "\n";
        }
    }

    destroy_environment(atm, bas, env);
    destroy_basis(basis);
    destroy_molecule(geom);

    end_tests(h);

    return h.exitcode;
}